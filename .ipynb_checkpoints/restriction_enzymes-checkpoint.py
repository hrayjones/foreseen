from Bio.Seq import Seq


class RestrictionEnzyme(Seq):
    def __init__(self,
                 enzyme_name:str,
                 enzyme_seq:str,
                 head_after_digestion:str,
                 tail_after_digestion:str
                 ):
        super().__init__(enzyme_seq)
        self.enz_name = enzyme_name
        self.head_add = - len(tail_after_digestion)
        self.tail_add = len(head_after_digestion)
        self.length = len(enzyme_seq)

def csv_to_enzyme_list(        
        csv_file: str,
):
    import pandas as pd
    from Bio.Seq import Seq
    
    ### read csv 
    df = pd.read_csv(csv_file, sep=';')
    df
    ### First case with enzymes that do not have spaces between the recognized sequence: 
    # 1. first remove the Enzymes that recognize trough spaces / remove substring N  
    # filter the rows that contain the substring
    substring = 'N'
    filter = df['Sequence'].str.contains(substring)
    filtered_df = df[~filter]
    filtered_df
    ### We need now to create 2 columns for the enzumes, the sequence without ( ^ v) and the pieces that are kept, the tail and head. 
    ### for the tail and head we need to take into account only the ^ because thats the one that indicates where the cut is done in the same strand
    ### that is recognized by the enzyme. 
        # 1. remove the v in another column 
    filtered_df['Sequence_wo_cuts'] = filtered_df['Sequence'].str.replace('v', '')
        #print(filtered_df)
        # 2. we divide the sequence then by the ^ Use the str.split method to split the data based on ^
    filtered_df[['pt1', 'pt2']] = filtered_df['Sequence_wo_cuts'].str.split('^', expand=True) #part 1 and part 2. 
        #print(filtered_df)
        # 3. Finish the removal of non ATCG characters 
    filtered_df['Sequence_wo_cuts'] = filtered_df['Sequence_wo_cuts'].str.replace('^', '')
    #print(filtered_df)
        # 4. Now depending on the Extreme from which the enzyme reads: 5' or 3' we define pt1 or pt2 as head or tail for the enzyme. The case for blunt end 
        #    is the same as 5'. For now, to have a first working example we will use only the 5' and blunt. 
    filtered_df = filtered_df[~filtered_df["Extreme"].str.contains("3", na=False)]
    #print(filtered_df)
        #
    filtered_df = filtered_df.reset_index()  # make sure indexes pair with number of rows
    ######## 
    ### restriction enzyme class designed by Lunkyadi Sucipto, the class includes an enzyme name, sequence, 
    ### head that is pt1 and tails that is pt2. 

    ### list of enzymes: 
    enzyme_list = []
    #
    for index, row in filtered_df.iterrows():
        # a. create temporal 
        tmp = RestrictionEnzyme(enzyme_name=row['Name'],
                                enzyme_seq=row['Sequence_wo_cuts'],
                                head_after_digestion=row['pt1'],
                                tail_after_digestion=row['pt2'])
        # b. look at what is added 
        #print(row['Name'], row['Sequence_wo_cuts'],row['pt1'], row['pt2'])
        # c. add to list 
        enzyme_list.append(tmp)
    
    #print(len(enzyme_list))
    #print(enzyme_list[1].enz_name,enzyme_list[1].head_add,enzyme_list[1].tail_add,enzyme_list[1].length)    
    # result: list of enzymes 
    return enzyme_list
