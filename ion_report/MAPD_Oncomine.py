__author__ = "Jonathan Serrano"
__version__ = "1.5.0"
__date__ = "11-15-2022"

import argparse
import json
import os
import subprocess
import warnings
import pandas as pd
import requests
import urllib3
import configparser

warnings.filterwarnings('ignore', category=UserWarning, module='openpyxl')
warnings.simplefilter(action='ignore', category=FutureWarning)
urllib3.disable_warnings()

def get_options():
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", "--config", type=str, required=True,
        help="Path to config file (e.g. ion_config.conf)")
    parser.add_argument("-r", "--runid", type=str, required=True,
                        help="Run ID 22-MGON37")
    parser.add_argument("-i", "--input", type=str, required=False, default="/mnt/Z_drive/Molecular/IonTorrent/oncosolid_autoreport/dropoff/",
                        help="Template Sheet Input Directory")
    parser.add_argument("-o", "--output", type=str, required=False, default="/home/ionadmin/ion_report/",
                        help="Output directory")
    return parser.parse_args()

def load_config(config_path):
    config = configparser.ConfigParser()
    config.read(config_path)
    HOST=config["DEFAULT"]["HOST"]
    TOKEN=config["DEFAULT"]["TOKEN"]
    return HOST, TOKEN


def get_mapd_values(sampleName):
    mapd_values = {}
    mapd_values.setdefault("name", [])
    mapd_values.setdefault("mapd", [])
    print("Getting Sample: " + sampleName)
    apiUrl = f'https://{HOST}/api/v1/qcreport?fomat=json'
    print(apiUrl)
    authToken = TOKEN
    try:
        resp = requests.get(
            apiUrl,
            headers={'Content-Type': 'application/x-www-form-urlencoded',
                     'Authorization': authToken},
            params={'name': sampleName}, verify=False)
        r = json.loads(resp.text)
        qc_metrics = r[0]["qc_metrics"]
        for sam in qc_metrics:
            mapdVal = qc_metrics[sam]["MAPD"]
            mapd_values["name"].append(sam)
            mapd_values["mapd"].append(mapdVal)
    except:
        print("Could not API Call sample: " + sampleName)

    return (mapd_values)


def get_sample_list(template_file):
    ts = pd.read_excel(
        template_file, 'DNA', skiprows=5,
        index_col=None, na_values=['NA'], engine='openpyxl')
    barcode = ts['Bar code']
    lastNrow = barcode.count() - 1
    print("~~~~~~Last row read in worksheet: ", lastNrow)
    workSheet = ts.loc[:lastNrow, :].copy()
    workSheet['Sample_ID'] = workSheet['Accession #'].astype(str)
    return (workSheet['Accession #'])
def float_covert(x):
    try:
        return float(x)
    except ValueError:
        return None

def build_mapd_csv(sample_sheet):
    mapd_df = pd.DataFrame()
    for sam in sample_sheet:
        jj = get_mapd_values(sam)
        df = pd.DataFrame.from_dict(jj, orient='index').transpose()
        df = df.fillna(method='ffill')
        mapd_df = pd.concat([mapd_df, df], ignore_index=True)
    mapd_df = mapd_df[mapd_df["name"].str.contains("RNA") == False].copy()
    #print(mapd_df.to_string())
    #cutoff = mapd_df["mapd"].apply(lambda x: float(x))
    cutoff = mapd_df["mapd"].apply(float_covert)
    larger = cutoff > 0.50
    mapd_df["Flagged"] = larger
    print(mapd_df.to_string())
    return (mapd_df)


# def download_rscript():
#     ''' Function downloads and chmod permissions of the Rscript '''
#     linkname = "https://raw.githubusercontent.com/NYU-Molecular-Pathology/oncomine"
#     filename = "MGON-gene_coverage_plots.R"
#     cmd1 = ["curl -# -L", linkname, ">" + os.getcwd() + "/" + filename]
#     subprocess.call(cmd1)
#     cmd2 = ["chmod +rwx", os.getcwd() + "/" + filename]
#     subprocess.call(cmd2)


# def execute_rscript(args):
#     '''Function executes the local Rscript passing the args Rscript -vanilla PATH runID OUTPUT '''
#     filename = "MGON-gene_coverage_plots.R"
#     command = ["Rscript", "--vanilla", os.getcwd() + "/" + filename,
#                args.runid]
#     subprocess.call(command)


def main():
    args = get_options()
    HOST, TOKEN = load_config(args.config)
    xlsxFi = os.path.join(args.input, args.runid + '.xlsm')
    print("Reading file: " + xlsxFi)
    template_file = pd.ExcelFile(xlsxFi, engine='openpyxl')
    sample_sheet = get_sample_list(template_file)
    mapd_df = build_mapd_csv(sample_sheet)
    outFile = args.runid + "_mapd_values.csv"
    outDirPath = os.path.join(args.output, outFile)
    print("File output: " + outDirPath)
    mapd_df.to_csv(outDirPath, index=False)
    # download_rscript()
    # execute_rscript(args)


if __name__ == "__main__":
    main()
