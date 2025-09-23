#%%
import csv
import warnings
import matplotlib.pyplot as plt

class Patient: 

    all_patients = []

    death_age = []

    APOEgeno = {}

    def __init__(self, DonorID, APOEgeno: str, age_symp_on: float, age_diag: float, sex: str):

        self.DonorID = DonorID
        self.ABeta40 = None
        self.ABeta42 = None
        self.tTau = None
        self.pTau = None
        self.sex = sex
        self.death_age = None
        self.ed_lvl = None
        self.cog_stat = None
        self.age_symp_on = age_symp_on
        self.age_diag = age_diag
        self.head_inj = None
        self.thal_score = None
        self.APOEgeno = APOEgeno
        Patient.all_patients.append(self)

    def __repr__(self):
        return f"{self.DonorID} | sex: {self.sex} | APOEgeno {self.APOEgeno} | age_symp_on {self.age_symp_on} | age_diag {self.age_diag} | Death Age {self.death_age}"

    def get_id(self):
        return self.DonorID

    def get_APOEgeno(self):
        return self.APOEgeno
    
    def get_age(self):
        return self.age_symp_on
    def get_age_diag(self):
        return self.get_age_diag
    
    def get_death_age(self):
        return self.death_age


    @classmethod
    def combine_data(cls, filename: str):
            with open(filename, encoding="utf8") as f:
                reader = csv.DictReader(f)
                rows_of_patients = list(reader)
                #for line in csv create object
                for row in range(len(rows_of_patients)):
                    if Patient.all_patients[row].DonorID == rows_of_patients[row]["Donor ID"]:
                        if rows_of_patients[row]["Sex"] != "":
                            Patient.all_patients[row].sex = rows_of_patients[row]["Sex"]

                        if rows_of_patients[row]["Age at Death"] != "":
                            Patient.all_patients[row].death_age = int(rows_of_patients[row]["Age at Death"])

                        if rows_of_patients[row]["Highest level of education"] != "":
                            Patient.all_patients[row].ed_lvl = rows_of_patients[row]["Highest level of education"]

                        if rows_of_patients[row]["Cognitive Status"] != "":
                            Patient.all_patients[row].cog_stat = rows_of_patients[row]["Cognitive Status"]

                        if rows_of_patients[row]["Age of onset cognitive symptoms"] != "":
                            Patient.all_patients[row].age_symp_on = int(rows_of_patients[row]["Age of onset cognitive symptoms"])

                        if rows_of_patients[row]["Age of Dementia diagnosis"] != "":
                            Patient.all_patients[row].age_diag = int(rows_of_patients[row]["Age of Dementia diagnosis"])

                        if rows_of_patients[row]["Known head injury"] != "":
                            Patient.all_patients[row].head_inj = rows_of_patients[row]["Known head injury"]

                        if rows_of_patients[row]["Thal"] != "":
                            Patient.all_patients[row].thal_score = int(rows_of_patients[row]["Thal"].split()[1])
            
                    else:
                        warnings.warn("IDs do not match.")
   
   #     @classmethod
    #def instantiate_from_csv(cls, filename: str, other_file: str):
    #open csv and create list of all rows
        # with open(filename, encoding="utf8") as f:
        #     reader = csv.DictReader(f)
        #     rows_of_patients = list(reader)
        #     #for line in csv create object
        #     for row in rows_of_patients:
        #         Patient(
        #             DonorID = row['Donor ID'],
        #             APOEgeno = float(row['ABeta40 pg/ug']),
        #             ABeta42 = float(row['ABeta42 pg/ug']),
        #             tTau = float(row['tTAU pg/ug']),
        #             pTau = float(row['pTAU pg/ug'])
        #         )
        #     Patient.all_patients.sort(key = Patient.get_id)
        #     Patient.combine_data(other_file) 

#%% 
#from patient import Patient

#Patient.instantiate_from_csv('UpdatedLuminex.csv', 'UpdatedMetaData (2).csv')

for patient in Patient.all_patients:
    print(patient)
    
# %%
