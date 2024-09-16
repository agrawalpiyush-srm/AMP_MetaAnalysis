### Python code ####
#### Boxplot GC distribution #####

df = pd.read_csv('boxplot_input_file.csv')
order = ['PathExt','DEGs']

sns.set_style('whitegrid')
ax= sns.boxplot(data=df, order=order, notch=True)
ax= sns.stripplot(data=df, order=order)
ax.set(ylim=(0,1))

add_stat_annotation(ax, data=df, order=order, box_pairs=[("Up", "Down")], test='Mann-Whitney', text_format='full', loc='outside', verbose=2)

ax.set_ylabel("3'UTR GC Content(%)",fontsize=15,weight='bold')
ax.set_xlabel("Cancer Type",fontsize=15,weight='bold')
ax.figure.set_size_inches(5,5)

plt.savefig('Boxplot_Ouput.png', dpi=600)
