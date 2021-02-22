


**Input1: Expression_matrix:**  
A tab separated matrix, rows = features, columns = samples

|  | sample1 | sample2 | sample3 |
| ------------- | ------------- | ------------- | ------------- |
| featureA  | 5.5  | 10  | 3  |
| featureB  | 0.12  | 0  | 4  |

**Input2: Annotation_Table:**  
A tab separated table, with the same number of rows as columns in the Expression matrix and a column containing the categories(classes) to perform your differential expression analysis.
|  | some_annotation| your_col_label | another_annotation |
| ------------- | ------------- | ------------- | ------------- |
| sample1  | 5.5  | classA  | 3  |
| sample2  | 0.12  | classB  | 4  |
| sample3  | 0.12  | classA  | 4  |

**Input3: code_path:**

The path were the code where your script (and its libraries) is.

**Input 4** is the path were your want to save the results.

**Input 5** the name of the column in your Annotation_Table that contain the information of your categories. 

**Input6** the log Fold change (base 2) to make your treshold

**Input7** your cut-off adjusted p.value

**Input8** If you want a label for your results.
