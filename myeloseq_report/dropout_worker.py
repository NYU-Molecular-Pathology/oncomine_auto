"""dropout_worker.py: Automated workflow for oncomine myeloid assay."""
__author__      = "Kelsey Zhu"
__copyright__   = "Copyright 2022, Langone Pathlab"

import pandas as pd
import os
import configparser
import glob
import re

class dropout(object):
    def __init__(self, conf_file):
        config = configparser.ConfigParser()
        config.read(conf_file)
        self.WORK_DIR = config['MYELOSEQ']['MYELOSEQ_HOME']
        self.DROPOUT_DIR = config['MYELOSEQ']['DROPOUT_DIR']
        self.TARGET_BED = config['MYELOSEQ']['TARGET_BED']
        self.REFSEQ = config['MYELOSEQ']['REFSEQ']
        self.TRANS_LOOKUP = config['MYELOSEQ']['TRANS_LOOKUP']
        self.STYLE_CSS = config['MYELOSEQ']['STYLE_CSS']
        # read data files
        self.target_bed = pd.read_csv(self.TARGET_BED, sep="\t", skiprows=1)
        self.target_bed.drop_duplicates(subset=['target'], inplace=True)
        self.refseq_df = pd.read_csv(self.REFSEQ, sep="\t", skiprows=0)
        self.myeloseq_trans = pd.read_csv(self.TRANS_LOOKUP, sep="\t", skiprows=0)
        self.myeloseq_trans[["Transcript", "Suffix"]] = self.myeloseq_trans["Transcript"].str.split(".", n=1, expand=True)
        self.myeloseq_refseq = self.refseq_df.merge(self.myeloseq_trans, left_on='name', right_on='Transcript', how='inner')
        self._workbook = None

    @property
    def workbook(self):
        return self._workbook

    @workbook.setter
    def workbook(self, value):
        self._workbook = value

    def get_chip_samples(self, sample_sheet, chip=1):
        sample_sheet = pd.read_excel(sample_sheet, engine='openpyxl', sheet_name='DNA', skiprows=5)
        sample_sheet = sample_sheet[['Bar code', 'Accession #', 'DNA #', 'Chip #']]
        sample_sheet.dropna(subset=['Bar code', 'Chip #'], inplace=True)
        sample_sheet = sample_sheet[sample_sheet['Accession #'].notna()]
        sample_sheet['Chip #'] = sample_sheet['Chip #'].astype('int64')
        sample_sheet = sample_sheet[sample_sheet['Chip #'] == int(chip)]
        sample_sheet['Bar code'] = sample_sheet['Bar code'].astype('int64')
        sample_sheet['sample_id'] = sample_sheet['Accession #'] + "-" + sample_sheet['DNA #']
        sample_sheet['Barcode'] = sample_sheet.apply(lambda x: "IonXpress_00%s"%x['Bar code'] if x['Bar code'] < 10
                            else "IonXpress_0%s"%x['Bar code'],axis=1)
        return sample_sheet

    def check_amplicons(self, row):
        for k,v in row.items():
            if int(v) < 500:
                return "FAIL"
        return "PASS"

    def get_NTC_colname(self, names):
        for n in names:
            if 'NTC' in n:
                return n
        return None

    def get_refseq_trans(self, row):
        for idx, mrow in self.myeloseq_refseq.iterrows():
            if row['Gene'] != mrow['name2'] or row['end'] < mrow['cdsStart'] or row['start'] > mrow['cdsEnd']: continue
            if row['start'] >= mrow['cdsStart'] and row['end'] <= mrow['cdsEnd'] or \
                row['end'] >= mrow['cdsStart'] and row['start'] <= mrow['cdsStart'] or \
                    row['start'] <= mrow['cdsEnd'] and row['end'] >= mrow['cdsStart']:
                return mrow['name']
            else:
                return None

    # def remove_zeros(self, df):
    #     drop_cols = df.columns[(df == 0).sum() > 0.25 * df.shape[1]]
    #     #drop_cols = df.columns[(df < 500).sum() > 30]
    #     print(df.columns[(df.columns < 500).sum() < 30])
    #     df.drop(drop_cols, axis=1, inplace=True)
    #     return df
    def remove_zeros(self, df):
        numeric_df = df.select_dtypes(include=['number'])
        # Columns where more than 25% of values are zeros
        drop_cols_zero = numeric_df.columns[(numeric_df == 0).sum() > 0.25 * numeric_df.shape[0]]
        # Drop columns where the number of values less than 500 is more than 30
        drop_cols_lt_500 = numeric_df.columns[(numeric_df < 500).sum() > 35]
        all_drop_cols = drop_cols_zero.union(drop_cols_lt_500)
        dropped_samples = all_drop_cols.tolist()
        df.drop(all_drop_cols, axis=1, inplace=True)
        
        return df, dropped_samples



    def make_clickable(self, row):
        if row['Transcript_ID'] == None:
            name =  row['Gene']
        else:
            name = row['Transcript_ID']
        return '<a target="_blank" href="https://genome.ucsc.edu/cgi-bin/hgTracks?org=human&position={0}&db=hg19">{1}</a>'.format(
                        "%s:%s-%s"%(row['chrom'],row['start'],row['end']), name)

    def get_run_id(self):
        header = pd.read_excel(self.workbook, engine='openpyxl', sheet_name='DNA', nrows=1, skiprows=2)
        return list(header.columns.values)[2].replace("-DNA", "")

    def process_dropout(self, df_sheet_all, chip=1):
        df_sheet_all, dropped_samples = self.remove_zeros(df_sheet_all)
        #df_sheet_all = self.remove_zeros(df_sheet_all)
        sample_sheet = self.get_chip_samples(self.workbook, chip)
        matched_dropped_samples_ids = sample_sheet[sample_sheet['Barcode'].isin(dropped_samples)]['sample_id'].tolist()
        df_sheet_all.sort_index(inplace=True, axis=1)
        old_col_names = list(df_sheet_all.columns.values)
        old_col_names.remove("Gene")
        old_col_names.remove("Target")
        df_sheet_all.set_index(['Gene', 'Target'], inplace=True)
        df_sheet_all.rename(
            columns={i:j for i,j in zip(sample_sheet['Barcode'],sample_sheet['sample_id'])}, inplace=True
        )
        df_sheet_all["filter"] = df_sheet_all.apply(lambda x: self.check_amplicons(x),axis=1)
        df_sheet_all = df_sheet_all.loc[df_sheet_all["filter"] == "FAIL"]
        #print(df_sheet_all)
        df_sheet_all.drop(['filter'], axis=1, inplace=True)
        df_sheet_all.reset_index(inplace=True)
        anno_amp_df = self.target_bed.merge(df_sheet_all,left_on='target', right_on='Target',how='inner')
        anno_amp_df.insert (4, "Hotspot", anno_amp_df.apply(lambda x: 'yes' if 'Hotspot' in x['info'] else 'no', axis=1))
        anno_amp_df.drop(['target', 'strand','info'], axis=1, inplace=True)
        self.myeloseq_refseq.sort_values(by=['cdsStart','cdsEnd'],inplace=True)
        anno_amp_df.sort_values(by=['start','end'],inplace=True)
        anno_amp_df.insert (5, "Transcript_ID", anno_amp_df.apply(lambda x: self.get_refseq_trans(x), axis=1))
        anno_amp_df.insert(6, "Transcript", anno_amp_df.apply(lambda x: self.make_clickable(x), axis=1))
        anno_amp_df.rename(columns={'chrom':'Chrom','start':'Start','end':'End'}, inplace=True)
        return anno_amp_df, matched_dropped_samples_ids

    def shift_row_to_top(self, df, row_to_shift):
        """Shift row, given by row_to_shift, to top of df."""

        idx = df.index.tolist()
        idx.pop(row_to_shift)
        df = df.reindex([row_to_shift] + idx)

        return df

    def start(self):
        HTML_STR = """
<html>
  <head><title>MyeloSeq Dropouts</title></head>
  <link rel="stylesheet" type="text/css" href="{style_css}"/>
  <body>
    <h3>Dropped Samples (more than 35 loci < 500 coverage):</h3>
    <p>{dropped_samples}</p>
    <h3>Dropouts Table:</h3>
    {table}
  </body>
</html>
"""#%self.STYLE_CSS
#         HTML_STR = """
# <html>
#   <head><title>MyeloSeq Dropouts</title></head>
#   <link rel="stylesheet" type="text/css" href=%s/>
#   <body>
#     {table}
#   </body>
# </html>
# """%self.STYLE_CSS
        run_id = self.get_run_id()
        dropout_files = glob.glob(os.path.join(self.DROPOUT_DIR, run_id, "Auto_%s*.bcmatrix.xls*"%run_id))
        anno_amp_df_list = list()
        dropped_samples_all = []
        for dropout_file in dropout_files:
            try:
                m = re.search(r'Auto_%s-(\d)(.*)'%run_id, dropout_file)
                chip = m.group(1)
            except:
                chip = 1
            dropout = pd.read_table(os.path.join(self.DROPOUT_DIR,run_id,dropout_file))
            anno_amp_df, dropped_samples = self.process_dropout(dropout, chip)
            dropped_samples_all.extend(dropped_samples)
            #anno_amp_df = self.process_dropout(dropout, chip)
            anno_amp_df.set_index(['Chrom', 'Start','End','Hotspot','Gene','Transcript',
                                   'Transcript_ID','Target'], inplace=True)
            anno_amp_df = anno_amp_df.apply(lambda x: x)
            anno_amp_df_list.append(anno_amp_df)

        if len(anno_amp_df_list) > 1:
            result = anno_amp_df_list[0].join(anno_amp_df_list[1], how='outer')
        else:
            result = anno_amp_df_list[0]

        result[result > 500] = 0
        result = result.fillna(0)
        result.loc[('Total Number of Drop-outs', '0','0','no','NA','NA',
                                   'NA','NA')] = result.apply(lambda x: len(x[x > 1]))
        result = result.astype(int)
        result.reset_index(inplace=True)
        result = self.shift_row_to_top(result, result.shape[0] - 1)
        pd.set_option('colheader_justify', 'center')

        # OUTPUT AN HTML FILE
        res = result.copy()
        res = res.drop(columns=['Transcript_ID'])
        dropped_samples_str = ", ".join(dropped_samples_all)
        with open(os.path.join(self.WORK_DIR, "%s-dropouts.html" % run_id), 'w') as f:
            f.write(HTML_STR.format(style_css=self.STYLE_CSS, dropped_samples=dropped_samples_str, table=res.to_html(classes='mystyle', escape=False, index=False)))


        # with open(os.path.join(self.WORK_DIR, "%s-dropouts.html"%run_id), 'w') as f:
        #     f.write(HTML_STR.format(table=res.to_html(classes='mystyle',escape=False, index=False)))

        result.drop(columns=['Transcript'],inplace=True)
        result.to_csv(os.path.join(self.WORK_DIR,"%s-dropouts.tsv"%run_id), sep="\t", index=False)
