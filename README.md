# AMP_MetaAnalysis

In this study, we obtained the data for ~8.6 Million Antimicrobial Peptides (AMPs) from the AMPSphere database developed by Santos-Junior et.al.
We performed meta-analysis of these peptides to predict their class specific activities i.e., whether they are Antibacterial, Antifungal and Antivrial Peptides.

We used various in silico prediction tools to predict their activity and based on consensus prediction with high score, we classify them into specific class.
We also performed other analysis such as Amino acid composition, Physiochemical properties, Motifs Analysis, and Positional Preference of Residues.

Using this dataset, we also developed sevreal machine leraning models to predict their activity and we acheived very high AUROC.

Next, we further screened these peptides for various properties a potential therapeutic candidate should exhibit such as allergenicity, toxicity, cell penentration ability, half-life in blood, etc. and proposed list of candidates in each class with potential therapeutic properties.

Lastly, we selected top 10 peptides in each class, predicted their 3D structure and performed molecular docking studies using pathogenic protein from an organism representing each class.

# Machine Learning Models

We have provided the training, testing and valdiation data used to developed the Machine Learning Models. Here positive data consists of peptides, we screened after multiple in silico experiments, whereas negative dataset comprises of those peptides which have been used previously to develop class specific prediction methods.

The code, datasets and ML models for each class is provided in the MachineLearning folder. Instructions regarding their usage is provided in the README.md file present in the same folder
