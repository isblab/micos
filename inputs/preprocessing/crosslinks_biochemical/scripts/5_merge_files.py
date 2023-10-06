import pandas as pd
import sys

# df = pd.read_csv('mouse_BDP_xl.csv', header = None)
# df1 = pd.read_csv('mouse_PIR_xl.csv', header = None)
# df2 = pd.read_csv('mouse_DSSO_xl.csv', header = None)
# df3 = pd.concat([df,df1,df2], ignore_index = True)
# df3.to_csv('mouse_data.csv', index = False, header = None)
# #
# hd = pd.read_csv('human_BDP_xl.csv',header = None)
# hd1 = pd.read_csv('human_PIR_xl.csv',header = None)
# hd2 = pd.read_csv('human_DSSO_xl.csv',header = None)
# hd3 = pd.concat([hd,hd1,hd2], ignore_index = True)
# hd3.to_csv('human_data.csv', index = False, header = None)
# df3 = pd.DataFrame()
# for files in glob.glob('/home/muskaan/Documents/modeling_micos_complex/archive/choice_of_species/yeast/xlink/*'):
#     df = pd.read_csv(files)
#     df3 = pd.concat([df3, df], ignore_index = True)
#
# df3.to_csv('yeast_data.csv', index = False)


output_file = sys.argv[1]
concatenated_df = pd.DataFrame()

for input_file in sys.argv[2:]:
    df = pd.read_csv(input_file)
    concatenated_df = pd.concat([concatenated_df, df], ignore_index=True)

concatenated_df.to_csv(output_file, index=False)
