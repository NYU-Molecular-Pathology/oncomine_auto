"""myelo_worker.py: Automated workflow for oncomine myeloid assay."""
# Version 4.1 - compatible with new IR; BAM downloading module; processing only IR filtered VCF
__author__      = "Kelsey Zhu, Summer Yang"
__copyright__   = "Copyright 2022, Langone Pathlab"

import requests
import os
import glob
import pandas as pd
import numpy as np
import re
import ast
from time import time
import logging
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype
import configparser
import urllib3
import shutil
import sys
import zipfile
#from ion_automation.myeloidseq.dropout_worker import dropout
from dropout_worker import dropout

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("ion_reporter")
RESULTS = list()
class myeloseq(object):

    def __init__(self, conf_file):
        self.conf_file = conf_file
        config = configparser.ConfigParser()
        config.read(conf_file)
        self.HOST = config['DEFAULT']['HOST']
        self.TOKEN = config['DEFAULT']['TOKEN']
        self.UID = config['DEFAULT']['UID']
        self.MYELOSEQ_HOME = config['MYELOSEQ']['MYELOSEQ_HOME']
        self.VAR_HOME = config['MYELOSEQ']['VAR_HOME']
        self.WORK_DIR = config['MYELOSEQ']['WORK_DIR']
        self.DEST_PATH = config['MYELOSEQ']['DEST_PATH']
        self.AA_CODES = config['MYELOSEQ']['AA_CODES']
        self.MYELOSEQ_GENES = config['MYELOSEQ']['MYELOSEQ_GENES']
        self.INCL_FUNCS = config['MYELOSEQ']['INCL_FUNCS'].split(",")
        self.EXCL_CALLS = config['MYELOSEQ']['EXCL_CALLS'].split(",")
        self.LOCATIONS = config['MYELOSEQ']['LOCATIONS'].split(",")
        self.STYLE_CSS = config['MYELOSEQ']['STYLE_CSS']
        self.VAR_TYPES = config['SOLID']['VAR_TYPES'].split(";")
        self.BAM_DIR = config['DEFAULT']['BAM_DOWNLOADS_DIR']

        self.codon_df = pd.read_csv(self.AA_CODES, sep=',')
        self.myeloseq_68genes = pd.read_excel(self.MYELOSEQ_GENES, engine='openpyxl', sheet_name='68-gene', skiprows=0)
        self.nyu_myeloseq_50genes = pd.read_excel(self.MYELOSEQ_GENES, engine='openpyxl', sheet_name='50-gene', skiprows=0)

        self._workbook = None
        self._dropout = dropout(conf_file)

    @property
    def workbook(self):
        return self._workbook

    @workbook.setter
    def workbook(self, value):
        self._workbook = value

    def __del__(self):
        self._workbook = None

    def get_coverage(self, row):
        try:
            m = re.search(r'.*;DP=([0-9]+);.*', row['INFO'])
            return m.group(1)
        except:
            return None

    def get_AF(self, row):
        try:
            m = re.search(r'AF=(.+);AO=.*;TYPE=(.*);VARB=(.*);HS;.*', row['INFO'])
            af_list = m.group(1).split(",")
            if len(af_list) > 1:
                idx = af_list.index(max(af_list))
                return "%s:%s:%s:%s" %(row['ALT'].split(",")[idx], max(af_list),
                                       m.group(2).split(",")[idx], m.group(3).split(",")[idx])
            else:
                return "%s:%s:%s:%s" %(row['ALT'],m.group(1),m.group(2), m.group(3))
        except:
            return None

    def get_func_row(self, row):
        try:
            m = re.search(r'AF=.*;HS;FUNC=(.*)', row['INFO'])
            func = ast.literal_eval(m.group(1))[0]
            return "%s:%s:%s" %(func['transcript'], func['gene'],func['exon'])
        except:
            return None

    def get_tumor_AF(self, row):
        return row['allele_frequency_%']

    def is_artifact(self, row):
        if row['type'] == 'SNV':
            ref, alt = tuple(row['genotype'].split("/"))
            if ref != alt and row['ref'] not in row['genotype']:
                return True
        return False

    def get_ExAC_info(self, row):
        try:
            m = re.search(r'AMAF=(.+):GMAF=(.+):EMAF=(.+)', row['5000Exomes'])
            return "%s:%s:%s"%(m.group(1), m.group(2), m.group(3))
        except:
            return "NA:NA:NA"

    def get_read_depth(self, row):
        try:
            return row['allele_coverage'].split(",")[1].split("=")[1]
        except:
            return None

    def get_alt_maf(self, row):
        try:
            return min(row['maf'].split(":"))
        except:
            return row['maf']

    def is_hotspot(self, row):
        try:
            return "Hotspot" in row['Oncomine Variant Annotator v3.2']
        except:
            return False

    def is_mnv(self, row):
        try:
            return "[" in row['protein'] or "," in row['protein'] or ";" in row['protein']
        except:
            return False

    def get_download_link(self, sample):
        """
        Return a tuple (download_link, name) or None if not found.
        Iterates v1, v2, ... and returns the last valid one.
        """
        try:
            headers = {
                "Content-Type": "application/x-www-form-urlencoded",
                "Authorization": self.TOKEN,
            }
            i = 1
            last_ok = None

            while True:
                params = {"name": f"{sample}_v{i}", "format": "json"}
                r = requests.get(
                    f"https://{self.HOST}/api/v1/getvcf",
                    headers=headers,
                    params=params,
                    verify=False,
                    timeout=30,
                )
                if r.status_code == 200:
                    last_ok = r
                    i += 1
                else:
                    break

            if not last_ok:
                return None

            j = last_ok.json()

            # Pull out fields depending on structure
            if isinstance(j, list) and j:
                data_links = j[0].get("data_links")
                name_field = j[0].get("name")
            elif isinstance(j, dict):
                data_links = j.get("data_links")
                name_field = j.get("name")
            else:
                data_links = None
                name_field = None

            # Normalize data_links to a single string
            if isinstance(data_links, list):
                data_links = data_links[0] if data_links else None

            if data_links:
                return str(data_links), name_field
            else:
                return None
        except Exception:
            return None


    def download_zip(self, sample):
        """
        download the zip from the returned link, unzip locally.
        """
        try:
            print("downloading zip")
            download_link= self.get_download_link(sample)[0]

            if not download_link:
                return None, None, None
            download_link = download_link.replace(
            "https://DPZNKD3:443", f"https://{self.HOST}", 1)
            print("printing:", download_link)
            # Ensure destination directories exist

            # Decide a filename for the zip
            # Try Content-Disposition; else use URL; else fall back to sample.zip
            headers = {
                "Content-Type": "application/x-www-form-urlencoded",
                "Authorization": self.TOKEN}
            
            r = requests.get(download_link, headers=headers, timeout=120, allow_redirects=False, stream=True, verify=False)
            r

            os.makedirs(self.VAR_HOME, exist_ok=True) 
            output_zip=os.path.join(self.VAR_HOME, "temp.zip")

            with open(output_zip, "wb") as f:
                for chunk in r.iter_content(chunk_size=1024 * 1024):
                    if chunk:  # keep-alive chunks can be empty
                        f.write(chunk)
                print(f"Saved: {output_zip}")
            return output_zip

        except Exception as e:
            import traceback
            print("ERROR:", type(e).__name__, str(e), flush=True)
            traceback.print_exc()
            return None
    
    def get_tsv_file(self, sample):
        """
        From the unzipped directory, locate and return the three files (tsv, vcf, oncomine tsv).
        or (None, None, None) if not found
        """
        temp_dir = "temp"
        # os.makedirs(temp_dir, exist_ok=True)
        temp_zip = self.download_zip(sample)
        print(f"temp zip is {temp_zip}")
        # Step 1: unzip temp.zip into temp_dir
        with zipfile.ZipFile(temp_zip, 'r') as z:
            z.extractall(temp_dir)

        # Step 2: find nested zip file inside temp_dir
        nested_zip = None
        for root, _, files in os.walk(temp_dir):
            for f in files:
                if f.endswith(".zip"):
                    nested_zip = os.path.join(root, f)
                    break
            if nested_zip:
                print(f"Copying files {nested_zip} into {self.DEST_PATH}")
                shutil.copy(nested_zip, os.path.join(self.DEST_PATH, "downloads"))
                break
        if not nested_zip:
            raise FileNotFoundError("No nested zip found in temp.zip")

        # Step 3: unzip the nested zip
        nested_dir = os.path.join(temp_dir, "nested")
        os.makedirs(nested_dir, exist_ok=True)
        with zipfile.ZipFile(nested_zip, 'r') as z:
            z.extractall(nested_dir)

        # Step 4: find the *subdirectory* under "Variants"
        variants_dir = None
        target_subdir = None
        for root, dirs, _ in os.walk(nested_dir):
            if os.path.basename(root) == "Variants":
                if dirs:
                    # assume only one target directory inside Variants
                    target_subdir = os.path.join(root, dirs[0])
                    print(target_subdir)
                break
        if not target_subdir:
            raise FileNotFoundError("No subdirectory under Variants found")

        # Step 5: copy that subdir to downloads_dir
        target_dir = self.VAR_HOME
        shutil.copytree(
            target_subdir,
            os.path.join(target_dir, os.path.basename(target_subdir)))

        # Step 6: clean up
        os.remove(temp_zip)
        shutil.rmtree(temp_dir)

        print(f"Saved directory {target_subdir} to {target_dir}")

        sample_pair=os.path.basename(target_subdir)
        file_path = os.path.join(self.VAR_HOME,sample_pair)
        return glob.glob(os.path.join(file_path, "%s*-full.tsv" % sample_pair))[0],\
            glob.glob(os.path.join(file_path, "%s*_Filtered_*.vcf" %sample_pair))[0], \
            glob.glob(os.path.join(file_path, "%s*_Filtered_*-oncomine.tsv" % sample_pair))[0]
# adding BAM downloading 
    def download_bam_file(self, url:str, sample: str, run_id:str):
        """
        Download a BAM file from the given inputBam URL and save it as {sample_name}.bam.
        """
        headers = {
            "Content-Type": "application/x-www-form-urlencoded",
            "auth": self.TOKEN,
            "Connection": "close"
        }
        
        run_path = os.path.join(self.BAM_DIR, run_id)
        os.makedirs(run_path, exist_ok=True) 

        out_path = os.path.join(run_path, f"{sample}.bam")
        
        try:
            with requests.get(url, headers=headers, stream=True, verify=False, timeout=600) as r:
                r.raise_for_status()
                print(f"Downloading {sample} BAM files")
                with open(out_path, "wb") as f:
                    for chunk in r.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
            return out_path
        except requests.RequestException as e:
            print(f"Failed to download {sample}: {e}")
            return None

    def fetch_and_download_bams(self, sample: str, run_id: str):
        """
        Fetch inputBam links for an analysis and download the BAM files.
        
        Returns a dict mapping sampleName -> bam_path.
        """
        analysis= self.get_download_link(sample)[1]
        results = {}
        url = f"https://{self.HOST}/api/v1/getAssociatedBamfiles"
        params = {"name": analysis, "type": "analysis"}
        headers = {"Content-Type": "application/x-www-form-urlencoded", "Authorization": self.TOKEN}
        
        try:
            resp = requests.get(url, headers=headers, params=params, verify=False, timeout=120)
            resp.raise_for_status()
            data = resp.json()
            
            for item in data:
                for sample in item.get("sampleDetails", []):
                    if sample.get("sampleRole") != "dna":
                        continue

                    sample_name = sample.get("sampleName")
                    input_bams = sample.get("inputBam", [])
                    if input_bams:
                        bam_url = input_bams[0]  # usually one inputBam
                        bam_url = bam_url.replace(
            "https://DPZNKD3:443", f"https://{self.HOST}", 1)
                        bam_path = self.download_bam_file(bam_url, sample_name, run_id)
                        results[sample_name] = bam_path
            return results
        
        except requests.RequestException as e:
            print(f"Error fetching input BAMs: {e}")
            return {}

    def clean_up(self):
        rm_downloads_cmd = 'rm -f %s/*.zip' % os.path.join(self.DEST_PATH, "downloads")
        os.system(rm_downloads_cmd)

    def copy_files(self, run_id):
        logger.info("Generating QC plots...")
        QC_plot_cmd = "Rscript MyeloSeq-QC-plot.R " + self.MYELOSEQ_HOME
        os.system(QC_plot_cmd)
        chrX_plot_cmd = "Rscript development_script_chrX_coverage.R " + run_id + " " + self.conf_file
        os.system(chrX_plot_cmd)

        logger.info("Copying the QC plots over to the Z drive...")
        mkdir_cmd = 'mkdir -p %s' % os.path.join(self.DEST_PATH, "reports/%s" % run_id)
        os.system(mkdir_cmd)

        plot_cp_cmd = 'cp -f %s %s' % ("*.pdf", os.path.join(self.DEST_PATH, "reports/%s" % run_id))
        os.system(plot_cp_cmd)

        csv_cp_cmd = 'cp -f %s %s' % ("*.csv", os.path.join(self.DEST_PATH, "reports/%s" % run_id))
        os.system(csv_cp_cmd)

        logger.info("Copying the report over to the Z drive")
        cp_cmd = 'cp -f %s %s' % ("%s.xlsx" % run_id, os.path.join(self.DEST_PATH, "reports/%s" % run_id))
        os.system(cp_cmd)

        if os.path.exists("%s-dropouts.html" % run_id):
            html_cmd = 'cp -f %s %s' % ("%s-dropouts.html" % run_id,
                                        os.path.join(self.DEST_PATH, "reports/%s" % run_id))
            os.system(html_cmd)

        sc_cp_cmd = 'cp -f %s %s' % ("sc_filtered_variants.tsv",
                                    os.path.join(self.DEST_PATH, "reports/%s/%s_SC_Variants.tsv" %
                                                (run_id, run_id)))
        os.system(sc_cp_cmd)

        if os.path.exists("sc2_filtered_variants.tsv"):
            sc2_cp_cmd = 'cp -f %s %s' % ("sc2_filtered_variants.tsv",
                                        os.path.join(self.DEST_PATH, "reports/%s/%s_SC2_Variants.tsv" %
                                                    (run_id, run_id)))
            os.system(sc2_cp_cmd)

        # copy variant tsv/vcf files over to the Z drive
        logger.info("Copying variant files over to the Z drive")
        if os.path.exists(os.path.join(self.DEST_PATH, "downloads/%s" % run_id)):
            rm_cmd = 'rm -rf %s' % os.path.join(self.DEST_PATH, "downloads/%s" % run_id)
            os.system(rm_cmd)

        cp_var_cmd = 'cp -r %s %s' % (self.VAR_HOME,
                                    os.path.join(self.DEST_PATH, "downloads/%s" % run_id))
        os.system(cp_var_cmd)

        rm_var_cmd = 'rm -rf %s/*' % self.VAR_HOME
        os.system(rm_var_cmd)

        rm_outfiles_cmd = 'rm -f %s/*.csv %s/*.pdf %s/*.xlsx %s/*.tsv %s/*.html' % (
            self.MYELOSEQ_HOME, self.MYELOSEQ_HOME, self.MYELOSEQ_HOME,
            self.MYELOSEQ_HOME, self.MYELOSEQ_HOME)
        os.system(rm_outfiles_cmd)

        self.clean_up()

    def write_to_excel(self, df):
        try:
            df = df.fillna("")
            run_id = df.iloc[0,0]
            logger.info("run ID: %s"%run_id)
            # Create a Pandas XlsxWriter engine.
            writer = pd.ExcelWriter("%s.xlsx"%run_id, engine='xlsxwriter')
            df.to_excel(writer, sheet_name='DNA', index=False)
            workbook = writer.book

            # Light red fill with dark red text.
            format11 = workbook.add_format({'bg_color': '#FFC7CE',
                                           'font_color': '#9C0006'})

            # Light yellow fill with dark yellow text.
            format22 = workbook.add_format({'bg_color': '#FFEB9C',
                                           'font_color': '#9C6500'})

            # Green fill with dark green text.
            format33 = workbook.add_format({'bg_color': '#C6EFCE',
                                           'font_color': '#006100'})

            worksheet = writer.sheets['DNA']

            # Add some cell formats.
            format1 = workbook.add_format({'num_format': '#,##0.00'})
            format2 = workbook.add_format({'num_format': '0%'})

            # Define the formats
            format6 = workbook.add_format({'bg_color': '#B9D3EE', 'border_color': '#BFBFBF', 'border': 1})  # dark blue
            format4 = workbook.add_format({'bg_color': '#FFFFFF', 'border_color': '#BFBFBF', 'border': 1})  # white
            format5 = workbook.add_format({'bg_color': '#E5F8F3', 'border_color': '#BFBFBF', 'border': 1})  # light green
            format3 = workbook.add_format({'bg_color': '#DBE8F0', 'border_color': '#BFBFBF', 'border': 1})  # light blue

            for row in range(1, df.shape[0] + 1):
                worksheet.write(row, 0, df.iloc[row-1, 0], format5)

            print (list(self.nyu_myeloseq_50genes['NYU.Myeloid']))
            print(list(self.myeloseq_68genes['Oncomine.Myeloid']))

            # Set the color for the starting row
            current_color = 'Dark'
            # Format the 1st row

            for column in range(1, df.shape[1]):  # format the first 2 columns
                worksheet.write(1, column, df.iloc[0, column], format3)

            # Start formatting from the 2nd row until the end of the df
            for row in range(2,df.shape[0]+1):
                # if the id of the row is the same as the id of the previous row
                if df.iloc[row - 1, 1] == df.iloc[row - 2, 1]:
                    if current_color == 'Dark':
                        sample_format = format3
                    elif current_color == 'Light':
                        sample_format = format4
                    for column in range(1, df.shape[1]):  # format the first 2 columns
                        worksheet.write(row, column, df.iloc[row - 1, column], sample_format)
                # if it's different than that of the previous row switch the colors
                else:
                    if current_color == 'Dark':
                        current_color = 'Light'
                    elif current_color == 'Light':
                        current_color = 'Dark'
                    for column in range(1, df.shape[1]):  # format the first 2 columns
                        worksheet.write(row, column, df.iloc[row - 1, column],
                                        format3 if current_color == 'Dark' else format4)
            # Set the column width and format.
            worksheet.set_column('B:B', 22, format1)
            worksheet.set_column('D:D', 18, format1)
            worksheet.set_column('N:N', 18, format1)
            worksheet.set_column('P:P', 18, format1)

            # Set the format but not the column width.
            worksheet.set_column('L:L', None, format2)

            for row in range(1, df.shape[0] + 1):
                if df.iloc[row-1, 5] != 'FUSION' and df.iloc[row-1, 4] != 'NA':
                    print(df.iloc[row - 1, 4], df.iloc[row - 1, 5])
                    if df.iloc[row-1, 4] in list(self.nyu_myeloseq_50genes['NYU.Myeloid']):
                        worksheet.write(row, 4, df.iloc[row - 1, 4], format33)
                    elif df.iloc[row-1, 4] in list(self.myeloseq_68genes['Oncomine.Myeloid']):
                        worksheet.write(row, 4, df.iloc[row-1, 4], format22)
                    else:
                        worksheet.write(row, 4, df.iloc[row-1, 4], format11)

            # Close the Pandas Excel writer and output the Excel file.
            writer.close()
        except:
            raise

    def get_vcf_fusion_key(self, row):
        try:
            m = re.search(r'(\w*)\[*\]*(chr\d+):(\d+)\[*\]*(\w)*',row['ALT'])
            return row['ID'].split("_")[0]
        except Exception as e:
            return None

    def get_tsv_fusion_key(self, row):
        try:
            if row['type'] == 'RNAExonVariant':
                if row['gene'] == 'EGFR|EGFR':
                    return 'EGFR-EGFR.E1E8.DelPositive.1'
                elif row['gene'] == 'MET|MET':
                    return 'MET-MET.M13M15'
            else:
                return row['# locus'].split("_")[1]
        except:
            return row['# locus']

    def get_fusion_read_counts(self, row):
        try:
            m = re.search(r'SVTYPE=([RNAExonVariant]*[Fusion]*);READ_COUNT=(.+);GENE_NAME=(.+);RPM=(.+);NORM_COUNT=(.+)',
                            row['INFO']);
            return m.group(2)
        except:
            return 0

    def get_fusion_RPM(self, row):
        try:
            m = re.search(r'SVTYPE=([RNAExonVariant]*[Fusion]*);READ_COUNT=(.+);GENE_NAME=(.+);RPM=(.+);NORM_COUNT=(.+)',
                          row['INFO']);
            return m.group(4)
        except:
            return 0

    def get_codon_letter(self, code):
        try:
            idx = self.codon_df.loc[self.codon_df['Codon'] == code].index
            return self.codon_df.iloc[idx]['Letter'].values[0]
        except:
            return None

    def get_codon_code(self, AA_change):
        try:
            return AA_change.replace("Ala", "A") \
                             .replace("Arg", "R") \
                             .replace("Asn", "N") \
                             .replace("Asp", "D") \
                             .replace("Cys", "C") \
                             .replace("Gln", "Q") \
                             .replace("Glu", "E") \
                             .replace("Gly", "G") \
                             .replace("His", "H") \
                             .replace("Ile", "I") \
                             .replace("Leu", "L") \
                             .replace("Lys", "K") \
                             .replace("Met", "M") \
                             .replace("Phe", "F") \
                             .replace("Pro", "P") \
                             .replace("Ser", "S") \
                             .replace("Thr", "T") \
                             .replace("Trp", "W") \
                             .replace("Tyr", "Y") \
                             .replace("Val", "V") \
                             .replace("Ter", "X")
        except:
            return AA_change

    def get_AA_Change(self, row):
        if isinstance(row['Amino Acid Change'], str):
            aa_change = row['Amino Acid Change'].split("|")[0]
            return self.get_codon_code(aa_change)
        else:
            return 'NA'

    def empty_row(self,row):
        try:
            if row['Type'] in ['SNV', 'INDEL'] and row['% Frequency'].astype(float) > 0:
                return False
            elif row['Type'] == 'RNAExonVariant':
                if row['Read Counts'] == 'NA' and row['Read/M'] == 'NA':
                    return True
                elif row['Genes'] in (['BRAF','EGFR']):
                    return True
                elif row['Read Counts'].isnull().values.any() or \
                    row['Read/M'].isnull().values.any():
                    return True
            elif row['Type'] in ['SNV', 'INDEL'] and str(row['Variant Effect']) == 'nan' \
                    and str(row['Transcript']) == 'nan' and str(row['Coding']) == 'nan':
                return True
            else:
                return False
        except:
            return False

    def select_gene(self, row):
        try:
            if row['type'] == 'FUSION':
                genes = row['gene'].split("|")
                exons = row['exon'].split("|")
                return "%s(%s)::%s(%s)" %(genes[0], exons[0], genes[1], exons[1])
            if row['gene'].endswith("|"):
                return row['gene'].replace("|","")
            elif not "|" in row['gene']:
                return row['gene']
            else:
                row_genes = row['gene'].strip().split("|")
                sel_gene = list(set(row_genes) & set(list(self.myeloseq_68genes['Oncomine.Myeloid'])))
                return sel_gene[0]
        except:
            return row['gene']

    def remove_function_bar(self, row):
        try:
            if row['function'].startswith("|") or row['function'].endswith("|"):
                return row['function'].replace("|", "")
            else:
                return row['function']
        except:
            return row['function']

    def get_location(self, row):
        try:
            return row['location'].split(":")[1]
        except:
            return row['location']
        
    def oncomine_in(self, row):
        try:
            value = row['Oncomine Variant Annotator v3.2']
            is_not_nan_and_has_length = pd.notna(value) and len(str(value)) > 0
            return is_not_nan_and_has_length
        except:
            return False

    def process_sample(self, args):
        sample, run_id, bar_code, logger = tuple(args)
        logger.info("start processing %s from %s", sample, run_id)

        filtered_tsv_path, filtered_vcf_path, oncomine_tsv_path = self.get_tsv_file(sample)
        print(f"Filtered_tsv_path is {filtered_tsv_path}")
        print(f"Filtered_vcf_path is {filtered_vcf_path}")

        if not filtered_tsv_path or not filtered_vcf_path:
            logger.warning("%s missing TSV/VCF paths", sample)
            return self._empty_result(sample, run_id, bar_code)

        ion_variants = None

        try:
            print("Reading the vcf FILTERED")
            filtered_vcf = pd.read_csv(filtered_vcf_path, sep="\t", skiprows=188)

            if not filtered_vcf.empty:
                filtered_vcf = filtered_vcf.loc[filtered_vcf['FILTER'].isin(["PASS", "."])]
                filtered_vcf["fusion_key"] = filtered_vcf.apply(self.get_vcf_fusion_key, axis=1)
                filtered_vcf["Read Counts"] = filtered_vcf.apply(self.get_fusion_read_counts, axis=1)
                filtered_vcf["Read/M"] = (
                    filtered_vcf.apply(self.get_fusion_RPM, axis=1)
                    .astype(float).astype(np.int64)
                )
                filtered_vcf["locus"] = (
                    filtered_vcf["#CHROM"].astype(str) + ":" + filtered_vcf["POS"].astype(str)
                )
                print(filtered_vcf.head(3).to_string())

                filtered_tsv = pd.read_csv(filtered_tsv_path, sep="\t", skiprows=2)
                filtered_tsv = filtered_tsv.loc[
                    filtered_tsv["filter"].isin(["PASS", "GAIN", ".", "LOSS"])
                    & (~filtered_tsv["type"].isin(self.EXCL_CALLS))
                ]
                pattern = "|".join(map(re.escape, filtered_vcf["locus"]))
                ion_variants = filtered_tsv[filtered_tsv["# locus"].str.contains(pattern)]

                filtered_vcf = filtered_vcf[["fusion_key", "Read Counts", "Read/M"]]

        except Exception as e:
            logger.error("Error processing %s: %s", sample, str(e))
            ion_variants = None

        # --- handle empty case early ---
        if ion_variants is None or ion_variants.empty:
            print("saving empty results")
            return self._empty_result(sample, run_id, bar_code)

        # Add annotations
        print("Adding annotations")
        ion_variants["tumor_AF"] = ion_variants.apply(self.get_tumor_AF, axis=1)
        ion_variants["ExAC_info"] = ion_variants.apply(self.get_ExAC_info, axis=1)
        ion_variants["DP"] = ion_variants.apply(self.get_read_depth, axis=1)
        ion_variants["MAF"] = ion_variants.apply(self.get_alt_maf, axis=1)
        ion_variants["HS"] = ion_variants.apply(lambda x: "yes" if self.is_hotspot(x) else "", axis=1)
        ion_variants["artifact"] = ion_variants.apply(self.is_artifact, axis=1)
        ion_variants["function"] = ion_variants.apply(self.remove_function_bar, axis=1)
        ion_variants["splice_site"] = ion_variants.apply(self.get_location, axis=1)
        ion_variants["MNV"] = ion_variants.apply(lambda x: "MNV" if self.is_mnv(x) else "", axis=1)
        ion_variants["fusion_key"] = ion_variants.apply(self.get_tsv_fusion_key, axis=1)
        ion_variants["gene"] = ion_variants.apply(self.select_gene, axis=1)
        ion_variants["# locus"] = ion_variants.apply(self.rename_locus, axis=1)
        ion_variants["anno_key"] = ion_variants.apply(
            lambda x: "%s:%s:%s" % (x["gene"], x["# locus"], x["protein"]), axis=1
        )

        new = ion_variants["ExAC_info"].str.split(":", n=2, expand=True)
        ion_variants.insert(0, "Run", run_id)
        ion_variants.insert(1, "Sample", sample)
        ion_variants.insert(2, "Barcode", bar_code)

        # print("Before filters are applied")
        # print(ion_variants.head(3).to_string())

        if not filtered_vcf.empty:
            ion_variants = ion_variants.merge(filtered_vcf, on="fusion_key", how="left")
        else:
            ion_variants["Read Counts"] = "NA"
            ion_variants["Read/M"] = "NA"
        print("Before filters are applied")
        print(ion_variants.head(3).to_string())
        # Drop unnecessary columns
        ion_variants.drop(
            columns=[
                "ExAC_info", "go", "5000Exomes", "hrun", "drugbank", "fusion_presence",
                "ratio_to_wild_type", "norm_count_within_gene", "filter",
                "allele_coverage", "allele_ratio", "pvalue", "dgv",
                "allele_frequency_%", "MyVariantDefaultDb_hg19", "phylop", "pfam",
                "location", "maf", "sift", "polyphen", "grantham", "normalizedAlt",
                "NamedVariants"
            ],
            inplace=True,
            errors="ignore"  # <-- ignore if missing
        )

        # Rename columns
        ion_variants.rename(
            columns={
                "# locus": "Locus", "gene": "Genes", "exon": "Exon",
                "transcript": "Transcript", "genotype": "Genotype",
                "type": "Type", "coding": "Coding", "function": "Variant Effect",
                "ref": "Ref", "tumor_AF": "% Frequency", "length": "Length",
                "protein": "Amino Acid Change", "coverage": "Coverage",
                "exac": "ExAC", "MAF": "ExAC_AF"
            },
            inplace=True
        )

        # Add new split columns
        ion_variants["AMAF"] = new[0]
        ion_variants["GMAF"] = new[1]
        ion_variants["EMAF"] = new[2]

        # Reorder columns
        ion_variants = ion_variants[
            ["Run", "Sample", "Barcode", "Locus", "Genes", "Type", "Exon",
            "Transcript", "Coding", "Variant Effect", "Genotype", "% Frequency",
            "ExAC_AF", "Amino Acid Change", "Read Counts", "Read/M",
            "AMAF", "GMAF", "EMAF", "HS", "Length", "Coverage",
            "Oncomine Variant Annotator v3.2"]
        ]
        
        logger.info("%s processed", sample)
        logger.info("After filters are applied")
        print(ion_variants.head(3).to_string())
        # Remove duplicates and bad rows
        # ion_variants = ion_variants.drop_duplicates(
        #     subset=["Sample", "Barcode", "Locus", "Genes", "Type", "Coding"]
        # )
        # ion_variants["bad"] = ion_variants.apply(self.empty_row, axis=1)
        # ion_variants = ion_variants.loc[ion_variants["bad"] != 1]


        return ion_variants


    def _empty_result(self, sample, run_id, bar_code):
        """Return a standard empty result DataFrame with NA values."""
        return pd.DataFrame({
            "Run": run_id, "Sample": sample, "Barcode": bar_code,
            "Locus": "negative", "Genes": "NA", "Type": "NA",
            "Exon": "NA", "Transcript": "NA", "Coding": "NA",
            "Variant Effect": "NA", "Genotype": "NA",
            "% Frequency": "NA", "ExAC_AF": "NA", "Amino Acid Change": "NA",
            "Read Counts": "NA", "Read/M": "NA", "AMAF": "NA",
            "GMAF": "NA", "EMAF": "NA", "HS": "NA", "Length": "NA",
            "Coverage": "NA", "Oncomine Variant Annotator v3.2": "NA"
        }, index=[0])


    def rename_genes(self, row):
        try:
            if row['Type'] == 'FUSION':
                genes = row['Genes'].split("|")
                exons = row['Exon'].split(",")
                return "%s(%s) - %s(%s)" %(genes[0],exons[0],genes[1],exons[1])
            else:
                if row['Genes'].endswith("|"):
                    return row['Genes'].replace("|", "")
                else:
                    return row['Genes'].replace("|", ",")
        except:
            return row['Genes']

    def rename_locus(self, row):
        try:
            if row['type'] == 'FUSION':
                if "_" in row['# locus']:
                    return row['# locus'].split("_")[0]
            return row['# locus']
        except Exception as e:
            logger.error(str(e))
            return row['# locus']

    def reform(self, colname,row): # reform exons to remove any "|"
        try:
            if row[colname].endswith("|"):
                return row[colname].replace("|", "")
            else:
                return row[colname].replace("|", ",")
        except:
            return row[colname]
    #clean up transcripts, amino acide code - only keeps the transcript (NM_xxxx), select correcponsing code and aa change
    def clean_transcript(self, transcript_str, coding_str, aa_str):
        if isinstance(transcript_str, str) and isinstance(coding_str, str) and isinstance(aa_str, str):
            transcript_splited = transcript_str.split("|")
            coding_splited = coding_str.split("|")
            aa_splited = aa_str.split("|")
            print(transcript_splited)
            print(coding_splited)
            nm_list = [s for s in transcript_splited if s.startswith('NM')]
            
            if not nm_list:  # If no strings start with 'NM'
                return None, None, None  # Return None for all the output
            
            selected_string = min(nm_list, key=len)  # Pick the shorter string - the shortest transcript naming
            selected_index = transcript_splited.index(selected_string)
            # Get the corresponding coding string based on the selected_index
            # If index is out of bounds for the coding_str, return None
            coding_selected = coding_splited[selected_index] if selected_index < len(coding_splited) else None
            aa_selected = aa_splited[selected_index] if selected_index < len(aa_splited) else None
            return selected_string, coding_selected, aa_selected
        else:
            return None, None, None

    def clean_exon(self, exon_str):
        if isinstance(exon_str, str):
            parts = exon_str.split("|")
            for part in parts:
                if part.isdigit():
                    return int(part)
        return None  

    def start(self):
        self.clean_up()
        RESULTS = list()
        try:
            ts = time()
            # read oncomine sample sheet
            logger.info(self.workbook)
            header = pd.read_excel(self.workbook, engine='openpyxl', sheet_name='DNA', nrows=1, skiprows=2)
            run_id = list(header.columns.values)[2].replace("-DNA","")
            sample_sheet = pd.read_excel(self.workbook, engine='openpyxl', sheet_name='DNA', skiprows=5)
            sample_sheet = sample_sheet.dropna(subset=[sample_sheet.columns[0]], how='all')
            sample_sheet['sample_id'] = sample_sheet['Accession #'] + "-" + sample_sheet['DNA #']
            logger.info(sample_sheet.to_string())  # shows headers with top 5 rows
            sc_sample_name = list(sample_sheet['sample_id'])[0]
            try:
                sc2_sample_name = list(filter(lambda x: "m2-Seraseq" in x, list(sample_sheet['sample_id'])))[0]
            except IndexError:
                sc2_sample_name = None
                print("No matches found for 'm2-Seraseq'")
            
            logger.info("SC sample: %s" %sc_sample_name)
            for sample, barcode in zip(list(sample_sheet['sample_id']), list(sample_sheet['Bar code'])):
                try:
                    if sample == "" or sample == None or str(sample) == 'nan': continue
                    df_result = self.process_sample([sample,run_id,barcode,logging.getLogger(sample)])
                    RESULTS.append(df_result)
                    logger.info("Processed sample %s -> %d rows", sample, len(df_result))
                    logger.debug(df_result.head(3).to_string()) 
                except:
                    pass
            if RESULTS and len(RESULTS) > 0:
                df = pd.concat(RESULTS)

                # Process SC
                sc_df = df.loc[df['Sample'] == sc_sample_name].copy()
                if not sc_df.empty:
                    sc_df.loc[:, 'Read Counts'] = sc_df.apply(
                        lambda x: 'NA' if x['Type'] in ['SNV', 'INDEL'] else x['Read Counts'], axis=1
                    )
                    sc_df.loc[:, 'Read/M'] = sc_df.apply(
                        lambda x: 'NA' if x['Type'] in ['SNV', 'INDEL'] else x['Read/M'], axis=1
                    )
                    sc_df.to_csv("sc_filtered_variants.tsv", index=False, sep="\t")

                # Process SC2
                if sc2_sample_name is not None:
                    sc2_df = df.loc[df['Sample'] == sc2_sample_name].copy()
                    if not sc2_df.empty:
                        sc2_df.loc[:, 'Read Counts'] = sc2_df.apply(
                            lambda x: 'NA' if x['Type'] in ['SNV', 'INDEL'] else x['Read Counts'], axis=1
                        )
                        sc2_df.loc[:, 'Read/M'] = sc2_df.apply(
                            lambda x: 'NA' if x['Type'] in ['SNV', 'INDEL'] else x['Read/M'], axis=1
                        )
                        sc2_df.to_csv("sc2_filtered_variants.tsv", index=False, sep="\t")

                # sample DF w/o SC1 SC2..
                sample_df = df.loc[~df['Sample'].isin([sc_sample_name, sc2_sample_name])].copy()
                sample_df.loc[:, 'Exon'] = sample_df['Exon'].apply(self.clean_exon)
                logger.info(sample_df.to_string())

                sample_df[['Transcript', 'Coding', 'Amino Acid Change']] = sample_df.apply(
                    lambda row: self.clean_transcript(row['Transcript'], row['Coding'], row['Amino Acid Change']),
                    axis=1, result_type='expand'
                )
                logger.info(sample_df.to_string())

                # Save to excel xlsx
                self.write_to_excel(sample_df)

                # alternate for final csv reporting
                add_df = sample_df.loc[~(sample_df['Sample'].str.contains("SC-DNA|NC-DNA", na=False))].copy()
                add_df = add_df.drop(columns=['Run','AMAF','GMAF','EMAF','Read Counts','Read/M','ExAC_AF'])
                add_df.loc[:, 'AA'] = add_df.apply(self.get_AA_Change, axis=1)

                add_df = add_df.rename(columns={
                    'Barcode': 'BC',
                    'Variant Effect': 'Variant.Effect',
                    'Amino Acid Change': 'Amino.Acid.Change',
                    '% Frequency': 'Frequency',
                    'HS': 'Info'
                })

                add_df = add_df[['Sample','BC','Locus','Genes','Type','Exon','Transcript','Coding','Variant.Effect',
                                'Genotype','Info','Length','Frequency','Amino.Acid.Change','AA','Coverage']]

                print(add_df.head(3).to_string())
                add_df.to_csv(f"{run_id}.csv", index=False, sep=",")

                # Continue downstream tasks
                self._dropout.workbook = self.workbook
                self._dropout.start()
                self.copy_files(run_id)

                for sample in list(sample_sheet['sample_id']):
                    try:
                        if sample == "" or sample == None or str(sample) == 'nan': continue
                        logger.info('Downloading %s sample BAM', sample)
                        self.fetch_and_download_bams(sample, run_id)
                    except:
                        pass
            logger.info('Took %s seconds to process samples', time() - ts)

        except Exception as e:
             logger.error(str(e))
