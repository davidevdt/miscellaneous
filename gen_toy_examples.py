### A small python program that allows generating 
### toy machine learning datasets. 

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from sklearn.datasets import make_moons
from sklearn.datasets import make_circles
import sys 



def plot_data(df, type_data=1, fig_size=(10,10), 
              mark_size=10, font_size=10): 
    
    ''' Function that takes a dataframe and the type of data 
    as input (along with other graphical parameters), and plots
    the corresponding data. 
    
    Parameters: 
        
    df (pandas.DataFrame): a dataframe, as the ones generated by 
                        make_data()
    
    type_data (int):    1 = classification, linear separation
                        2 = classification, linear non-separation
                        3 = classification, moons 
                        4 = classification, circles
                        5 = linear regression
                        6 = non-linear regression   
    
    fig_size (tuple): graphical parameter which defines the 
                        figure size. 
                        
    mark_size (int): graphical parameter whih defines the size 
                        of the plotted points. 
                        
    font_size (int): graphical parameter which defines the font 
                        size for the axes. 
                        
    
    Returns: 
        
    None
    
    '''
    
    
    if type_data not in (5, 6):  
        plt.figure(figsize=fig_size)
        plt.plot(df["x1"][df["y"]==0], df["x2"][df["y"]==0], "bx", markersize=mark_size)
        plt.plot(df["x1"][df["y"]==1], df["x2"][df["y"]==1], "ro", markersize=mark_size)
        plt.xlabel("x1", fontsize=font_size)
        plt.ylabel("x2", fontsize=font_size)
        plt.show()        

    else: 
        plt.figure(figsize=fig_size)
        plt.plot(df["x"], df["y"], "bo", markersize=mark_size)
        plt.xlabel("x", fontsize=font_size)
        plt.ylabel("y", fontsize=font_size)
        plt.show()





def export_data(df, filepath):
    
    '''Function that takes a dataframe and 
    a file path as an input, and exports the data
    frame into the path as a .csv file format. 
    
    Paramters: 
        
    df (pandas.DataFrame): a pandas dataframe. 
    
    filepath (str): path (including file name) for the  file that
    must be created. The file needs not include the '.csv' 
    extension.
    
    
    Returns: 
        
    None. 
    
    '''
    
    if not filepath[-3:] == 'csv': 
        filepath = filepath + '.csv' 
    
    df.to_csv(filepath, index=None, header=True)
           
    





def make_data(type_data=1, sample_size=20, starting_seed=2, 
                       noise_level = 0.1): 
    
    '''Function that takes as parameters the sample size,
    the desired level of noise, and the type of data that
    the user wishes to generate, and returns a dataset of 
    the generated data. The data are of the form X-->y, where
    X is two-dimensional in the classification case, and
    one-dimensional in the regression case. 
    
    Parameters: 
        
    type_data (int):    1 = classification, linear separation
                        2 = classification, linear non-separation
                        3 = classification, moons 
                        4 = classification, circles
                        5 = linear regression
                        6 = non-linear regression                    
    
    sample_size (size): number of samples to generate. When type_data = 1, 
                        the number of samples might be lower than the one 
                        specified depending on the noise_level parameter.
    
    starting_seed (int): seed to use for the random number generation
    
    noise_level (float): level of noise to use for type_data 3 and 4. 
                        Value of the standard deviation of the gaussian
                        when type_data = 5 and 6. When type_data = 1, 
                        the value specifies how much the points 
                        should lie far from the decision boundary. 
                        When type_data = 2, this value indicates how 
                        close to the boundary the values with flipped
                        labels should lie. 
                         
    
    Returns: 
        
    df (pd.DataFrame): the generated dataset. 
    
    
    '''
    
    
    if noise_level <= 0:
        print("Error: Noise Level should be greater than 0.")
        exit(1)
        

    if sample_size <= 0:
        print("Error: Sample Size should be greater than 0.")
        exit(1)
 

    if starting_seed <= 0:
        print("Error: Seed should be greater than 0.")
        exit(1)       
    
    
    #### Classification Models 
    if type_data in (1, 2): 
        # Generate predictors  
        np.random.seed( starting_seed )
        x1 = np.random.uniform(size=sample_size)
        x2 = np.random.uniform(size=sample_size)
    
        if type_data == 1: 
            # Case 1a: Linearly Separable Data, Perfect Separation 
            retain = (abs(x1-x2) > noise_level)
            x1_new = x1[retain]
            x2_new = x2[retain]
            new_size = sum(retain)
        
            y = np.zeros(new_size, dtype=int)
            y[x1_new < x2_new] = 1 
            # sum(y) --> Check 

            # Data Frame 
            df =  pd.DataFrame(dict(x1=x1_new, x2=x2_new, y=y))
    
    
        else: 
            # Case 1b: Lineary separable data, non-perfect separation  
            noise = (abs(x1-x2) < noise_level)
        
            y = np.zeros(sample_size, dtype=int)
            y[x1 < x2] = 1
            y[noise] = 1 - y[noise]
            # or : y[noise] = np.random.randint(2, size=sum(noise))
               
            # Data Frame
            df =  pd.DataFrame(dict(x1=x1, x2=x2, y=y))

    
    elif type_data in (3, 4):

        if type_data == 3:
            # Case 2a: non-linear case, moons
            X, y = make_moons(n_samples=sample_size, noise=noise_level, random_state=starting_seed)
        
            df = pd.DataFrame(dict(x=X[:,0], y=X[:,1], label=y))
            
            # Data Frame Columns 
            df.columns = ["x1", "x2", "y"]
            
        else:
            # Case 2b: non-linear case, circles
            X, y = make_circles(n_samples=sample_size, noise=noise_level, random_state=starting_seed)
            
            df = pd.DataFrame(dict(x=X[:,0], y=X[:,1], label=y))
                        
            # Data Frame Columns 
            df.columns = ["x1", "x2", "y"]

              
    #### Regression Models    
    else: 
        # Generate predictor
        np.random.seed( starting_seed )
        x1 = np.random.uniform(size=sample_size)
        
        if type_data == 5:
            # Case 1: Linear Regression
            b0 = 5.
            b1 = 2.
            
            y = b0 + b1*x1 + np.random.normal(scale=noise_level, size=sample_size)
            
            
            # Data Frame 
            df =  pd.DataFrame(dict(x=x1, y=y))

        else:
            x1 = np.random.uniform(0, 5, size=sample_size)
            
            b0 = 0.5
            b1 = 0.5
            b2 = 2.5
            b3 = -1
            
            y = b0 + b1*x1 + \
                b2*x1**2 + b3*x1**3 + \
                np.random.normal(scale=noise_level, size=sample_size)
            
     
            # Data Frame 
            df =  pd.DataFrame(dict(x=x1, y=y))

    
    return df 









if __name__ == "__main__":
    
    print("Welcome to the toy datasets generator program.")
    print("===============================================================")
    
    
    while True:
        selection_type = -1    
        
        while (int(selection_type) < 1 or int(selection_type) > 6 ):
            print("""Please select the type of dataset you want to generate,\nby choosing one of the following options:
                  
    [1] = classification, linear perfect separation
    [2] = classification, linear non-perfect separation
    [3] = classification, moons 
    [4] = classification, circles
    [5] = linear regression
    [6] = non-linear regression                 
    [9] = exit 
                  
              """)
                
            print("Select one of the options above:")
            read_input = input()
            selection_type = int( read_input )
            
            if selection_type not in (1, 2, 3, 4, 5, 6, 9):
                print("Please select options 1-6, or 9 to exit")    
            elif selection_type == 9:
                sys.exit()
            
            else:
                if selection_type == 1:
                    print("You selected 'classification, linear perfect separation'.")
                elif selection_type == 2: 
                    print("You selected 'classification, linear non-perfect separation'.")
                elif selection_type == 3: 
                    print("You selected 'classification, moons'.")
                elif selection_type == 4: 
                    print("You selected 'classification, circles'.")
                elif selection_type == 5: 
                    print("You selected 'linear regression'.")
                else: 
                    print("You selected 'non-linear regression'.")                
                # break
            
                
        sel_samp_size = -1
        sel_lev_noise = -1
        
        while sel_samp_size <= 0:
            print("Please select sample size:")
            read_input = input()
            sel_samp_size = int(read_input)
            
            if sel_samp_size <= 0 :
                print("Sample size must be larger than 0.")
            # else:
                # break
            
     
        while sel_lev_noise <= 0:
            if( selection_type not in (1, 5, 6)):
                print("Please select noise level:")
            elif selection_type == 1:
                print("Please select separation level:")
            else: 
                print("Please select error standard deviation:")
            read_input = input()
            sel_lev_noise = float(read_input)
            
            if sel_samp_size <= 0:
                print("Noise level must be larger than 0.")
            # else:
                # break
               
        
        data = make_data(type_data=selection_type, 
                         sample_size=sel_samp_size,
                         noise_level = sel_lev_noise)
        
        plot_data(data, type_data=selection_type, 
                  fig_size=(10,10), 
                  mark_size=10, 
                  font_size=10)
        
        print("""Data have been generated. Do you want to export it?
        [y] yes
        [q] exit without exporting data
        [r] restart 
              """)
        
        export_or_exit = input() 
        
        while export_or_exit not in ("y", "q", "r"):
            print("Please make your selection: y/q/r")
            export_or_exit = input()
            
        if export_or_exit == "q":
            print("Exiting the program...")
            sys.exit()
        elif export_or_exit == "y": 
            print("Please select the file destination, including filename:")
            filepath = input()
            print("Exporting the dataset...")        
            export_data(data, filepath)
            print("""Data exported. Exit or restart?
                [q] exit 
                [r] restart 
                     """)
                     
            restart_or_exit = input()
            
            while restart_or_exit not in ("q", "r"):
                print("Please make your selection: q/r")
                export_or_exit = input()
            
            if restart_or_exit == "q":
                sys.exit()
                
    
        print("===============================================================")   
        

    
