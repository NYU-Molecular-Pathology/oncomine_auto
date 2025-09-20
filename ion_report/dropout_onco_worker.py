"""dropout_worker.py: Automated workflow for oncomine focus assay."""
__author__      = "Kelsey Zhu"
__copyright__   = "Copyright 2022, Langone Pathlab"
__version__ = "1.0.2"

import pandas as pd
import os
import glob
import configparser
import re

class dropout(object):
    def __init__(self, conf_file):
        config = configparser.ConfigParser()
        config.read(conf_file)
        self.WORK_DIR = config['SOLID']['HOME_DIR']
        self.DROPOUT_DIR = config['SOLID']['DROPOUT_DIR']
        self.STYLE_CSS = config['DEFAULT']['STYLE_CSS']
        self._workbook = None

    @property
    def workbook(self):
        return self._workbook

    @workbook.setter
    def workbook(self, value):
        self._workbook = value

    def get_chip_samples(self, sample_sheet, chip=1):
        sample_sheet = pd.read_excel(sample_sheet, engine='openpyxl', sheet_name='DNA', skiprows=5)
        sample_sheet = sample_sheet[['Bar code', 'Accession #', 'DNA #', 'Chip # \n1 or 2']]
        sample_sheet.dropna(subset=['Bar code', 'Chip # \n1 or 2'], inplace=True)
        sample_sheet = sample_sheet[sample_sheet['Accession #'].notna()]
        sample_sheet['Chip #'] = sample_sheet['Chip # \n1 or 2'].astype('int64')
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

    def remove_zeros(self, df):
        drop_cols = df.columns[(df == 0).sum() > 0.25 * df.shape[1]]
        df.drop(drop_cols, axis=1, inplace=True)
        return df

    def get_run_id(self):
        header = pd.read_excel(self.workbook, engine='openpyxl', sheet_name='DNA', nrows=1, skiprows=2)
        return list(header.columns.values)[2].replace("-DNA", "")

    def process_dropout(self, df_sheet_all, chip):
        df_sheet_all = self.remove_zeros(df_sheet_all)
        sample_sheet = self.get_chip_samples(self.workbook, chip)
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
        df_sheet_all.drop(['filter'], axis=1, inplace=True)
        df_sheet_all.reset_index(inplace=True)
        return df_sheet_all

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
  <link rel="stylesheet" type="text/css" href=%s/>
  <body>
    {table}
  </body>
</html>
"""%self.STYLE_CSS
        run_id = self.get_run_id()
        dropout_files = glob.glob(os.path.join(self.DROPOUT_DIR,"Auto_%s*.bcmatrix.xls*"%run_id))
        anno_amp_df_list = list()
        for dropout_file in dropout_files:
            try:
                m = re.search(r'Auto_%s-(\d)(.*)'%run_id, dropout_file)
                chip = m.group(1)
            except:
                chip = 1
            #dropout_file = glob.glob(os.path.join(self.DROPOUT_DIR,"Auto_%s_*.bcmatrix.xls"%run_id))[0]
            dropout = pd.read_table(dropout_file)
            anno_amp_df = self.process_dropout(dropout, chip)
            anno_amp_df.set_index(['Gene', 'Target'], inplace=True)
            anno_amp_df = anno_amp_df.apply(lambda x: x)
            anno_amp_df_list.append(anno_amp_df)

        if len(anno_amp_df_list) == 2:
            result = anno_amp_df_list[0].join(anno_amp_df_list[1], how='outer')
        else:
            result = anno_amp_df_list[0]
        # anno_amp_df.set_index(['Gene','Target'], inplace=True)
        # result = anno_amp_df.apply(lambda x: x)
        result[result > 500] = 0
        result = result.fillna(0)
        result = result.astype(int)
        counts = result.apply(lambda x: len(x[x > 1]))
        result.reset_index(inplace=True)
        result.loc[result.shape[0]] = ['Total # Dropouts', ''] + counts.tolist()
        result = self.shift_row_to_top(result, result.shape[0] - 1)
        pd.set_option('colheader_justify', 'center')

        # OUTPUT AN HTML FILE
        res = result.copy()
        # change the html saved path to the dest dir maybe??
        with open(os.path.join(self.WORK_DIR, "%s-dropouts.html"%run_id), 'w') as f:
            f.write(HTML_STR.format(table=res.to_html(classes='mystyle',escape=False, index=False)))
        result.to_csv(os.path.join(self.WORK_DIR,"%s-dropouts.tsv"%run_id), sep="\t", index=False)

