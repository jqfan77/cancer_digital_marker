import tkinter as tk
import tkinter.filedialog as tkfd
import tkinter.ttk as ttk

import model as my_model

import os

window = tk.Tk()

window.title('Drug sensitivity prediction')
window.geometry("500x420")

###
L1 = 50
L2 = 250

U1 = 50
U3 = 125
U4 = 230
U5 = 350

### gene file
text_genefile_raw = '[No file]'

lbl_genefile = tk.Label(window, text='Gene data file: \n(gene symbols + FPKM/TPM) ', justify=tk.LEFT)
lbl_genefile.place(x=L1,y=U1)

lbl_genefile_2 = tk.Label(window, text='[No file]', fg='#A00000', justify=tk.LEFT)
lbl_genefile_2.place(x=L2,y=U1)


def clicked_genefile():
    file_genefile = tkfd.askopenfilename(initialdir=os.path.dirname(__file__))
    if len(file_genefile)>0:
        global text_genefile_raw
        text_genefile_raw = file_genefile
        if len(file_genefile)<30:
            lbl_genefile_2.config(text = file_genefile)
        else:
            lbl_genefile_2.config(text = '...'+file_genefile[-27:])
        lbl_genefile_2.config(fg = '#000000')
    else:
        text_genefile_raw = '[No file]'
        lbl_genefile_2.config(text = text_genefile_raw)
        lbl_genefile_2.config(fg = '#A00000')


btn_genefile = ttk.Button(window, text='Select file', command=clicked_genefile)
btn_genefile.place(x=L2, y=U1+20)

### tissue
lbl_tissue = tk.Label(window, text='Predict in tissue: ', justify=tk.LEFT)
lbl_tissue.place(x=L1, y=U3)

combo_tissue = ttk.Combobox(window)
combo_tissue['values'] = ('Colorectal')
combo_tissue.current(0)
combo_tissue.place(x=L2, y=U3)




'''### drug
lbl_drug = Label(window, text='Drug to predict: ')
lbl_drug.place(x=L1, y=U4)

combo_drug = Combobox(window)
combo_drug['values'] = ('1005_Cisplatin', '1073_5-Fluorouracil', '1080_Paclitaxel', '1512_Cyclophosphamide')
combo_drug_dict = {'1005_Cisplatin':1005, '1073_5-Fluorouracil':1073, '1080_Paclitaxel':1080, '1512_Cyclophosphamide':1512}
combo_drug.current(0)
combo_drug.place(x=L2, y=U4)'''

### drug display

drug_id_list = [1005, 1073, 1080, 1512]
drug_name_list = ['1005 Cisplatin', '1073 5-Fluorouracil', '1080 Paclitaxel', '1512 Cyclophosphamide']
lbl_drug_name_list = []
for i in range(len(drug_id_list)):
    lbl_drug_name = tk.Label(window, text=drug_name_list[i], justify=tk.LEFT)
    lbl_drug_name.place(x=L1, y=U4+i*20)
    lbl_drug_name_list.append(lbl_drug_name)

lbl_Y_predict_relative_title = tk.Label(window, text='z-score \n(relatively in \nthe tissue)', justify=tk.LEFT)
lbl_Y_predict_relative_title.place(x=L2-30, y=U4-54)
lbl_Y_predict_relative_list = []
for i in range(len(drug_id_list)):
    lbl_temp = tk.Label(window, text='???', justify=tk.LEFT)
    lbl_temp.place(x=L2-30, y=U4+i*20)
    lbl_Y_predict_relative_list.append(lbl_temp)

lbl_Y_predict_title = tk.Label(window, text='z-score \n(GDSC)', justify=tk.LEFT)
lbl_Y_predict_title.place(x=L2+65, y=U4-37)
lbl_Y_predict_list = []
for i in range(len(drug_id_list)):
    lbl_temp = tk.Label(window, text='???', justify=tk.LEFT)
    lbl_temp.place(x=L2+65, y=U4+i*20)
    lbl_Y_predict_list.append(lbl_temp)

lbl_Y_predict_log_ic50_title = tk.Label(window, text='log10(IC50/uM)', justify=tk.LEFT)
lbl_Y_predict_log_ic50_title.place(x=L2+130, y=U4-20)
lbl_Y_predict_log_ic50_list = []
for i in range(len(drug_id_list)):
    lbl_temp = tk.Label(window, text='???', justify=tk.LEFT)
    lbl_temp.place(x=L2+150, y=U4+i*20)
    lbl_Y_predict_log_ic50_list.append(lbl_temp)

    


### Run

lbl_result = tk.Label(window, text='Select the gene expression file, \nthen click `Predict`.', justify=tk.LEFT)
lbl_result.place(x=L2-80,y=U5)

def clicked_run():
    genefile_filename = text_genefile_raw
    
    '''select_dataset = combo_dataset_dict[combo_dataset.get()]'''
    select_dataset = 2 # no choice here

    tissue_selected_id_list = [14]  # no choice here

    global drug_id_list
    #drug_id_list = [1005, 1073, 1080, 1512]  # combo_drug_dict[combo_drug.get()]  # 改成1005,1073，1080, 1512 都要预测

    info_str, list_Y_predict_relative, list_Y_predict, list_Y_predict_log_ic50, list_color = my_model.predict(genefile_filename, select_dataset, tissue_selected_id_list, drug_id_list)

    for i in range(len(drug_id_list)):
        lbl_drug_name_list[i].config(fg = list_color[i])
        lbl_Y_predict_relative_list[i].config(text = list_Y_predict_relative[i])
        lbl_Y_predict_list         [i].config(text = list_Y_predict[i])
        lbl_Y_predict_log_ic50_list[i].config(text = list_Y_predict_log_ic50[i])

    lbl_result.config(text = info_str )

btn_run = ttk.Button(window, text='Predict', command=clicked_run)
btn_run.place(x=L1, y=U5)


### 

window.mainloop()

