import pandas as pd

class ImdbFile(file_name):
    def __int__(self):
        self.df = pd.read_table('imdb-reviews.tsv')

    def 




a = pd.read_table('au80.txt', header=None)
a = pd.read_table('au.txt', header=None)

s = [1, 0.901, 0.782]
den =(a.T/s).T.sum()


X = a.iloc[0]/den
X1 = np.log(X)
X2 = 1/X
X3 = 1/X**2
T = [10, 20, 40, 60, 80]



color_codes = ['#CD58E8', '#7F4CFF', '#FF5494', '#AC5745']
line_styles = ['-.', '--', ':', '-']
l =['CuSO4·5H2O|15K·min-1', 'CuSO4·5H2O|10K·min-1', 'CuSO4·5H2O|5K·min-1']


def filer(file_name):
    with open(file_name) as f:
        line_number = 0
        for line in f:
            line_number += 1
            if "Data=Data" in line:
                break
        file_data = pd.read_csv(file_name, skiprows=line_number, header=None)
        return file_data

data =[filer(f) for f in ['cu.txt', 'cu10.txt', 'cu5.txt',]]
fig, ax = plt.subplots()
for iteration in range(3):
    ax.plot(data[iteration][0],data[iteration][1], color = color_codes[iteration] ,linestyle=line_styles[iteration], label=l[iteration])
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.xaxis.tick_top()
ax.tick_params(direction='in')
plt.xlim(0,400)
pylab.legend(loc='lower right')












