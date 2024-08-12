import numpy as np
from numpy import linalg as LA
np.set_printoptions(threshold=np.inf)
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 15})
plt.rcParams.update({'text.usetex': True})
plt.rcParams.update({'text.latex.preamble' : r'\usepackage{amsmath}'})

t=1.0
#       1º region:
list_ferro_u1 = []
list_ferro_ne1 = []
with open('FERRO_U1.txt','r') as f:
    fu = [line[:-1] for line in f]

for i in range(len(fu)):
    list_ferro_u1.append(float(fu[i]))

with open('FERRO_ne1.txt','r') as f:
    fne = [line[:-1] for line in f]

for i in range(len(fu)):
    list_ferro_ne1.append(float(fne[i]))

ferro_tu1 = []
ferro_ne1 = []
for j in range(len(fu)):
    if (list_ferro_u1[j]!=-1 and list_ferro_ne1[j]!=-1):
        ferro_tu1.append(t/list_ferro_u1[j])
        ferro_ne1.append(list_ferro_ne1[j])
#       2º region:
list_ferro_u2 = []
list_ferro_ne2 = []
with open('FERRO_U2.txt','r') as f:
    fu = [line[:-1] for line in f]

for i in range(len(fu)):
    list_ferro_u2.append(float(fu[i]))

with open('FERRO_ne2.txt','r') as f:
    fne = [line[:-1] for line in f]

for i in range(len(fu)):
    list_ferro_ne2.append(float(fne[i]))

ferro_tu2 = []
ferro_ne2 = []
for j in range(len(fu)):
    if (list_ferro_u2[j]!=-1 and list_ferro_ne2[j]!=-1):
        ferro_tu2.append(t/list_ferro_u2[j])
        ferro_ne2.append(list_ferro_ne2[j])
########
#       1º region:
list_anti_u1 = []
list_anti_ne1 = []
with open('ANTI_U1.txt','r') as f:
    fu = [line[:-1] for line in f]

for i in range(len(fu)):
    list_anti_u1.append(float(fu[i]))

with open('ANTI_ne1.txt','r') as f:
    fne = [line[:-1] for line in f]

for i in range(len(fu)):
    list_anti_ne1.append(float(fne[i]))

anti_tu1 = []
anti_ne1 = []
for j in range(len(fu)):
    if (list_anti_u1[j]!=-1 and list_anti_ne1[j]!=-1):
        anti_tu1.append(t/list_anti_u1[j])
        anti_ne1.append(list_anti_ne1[j])
#       2º region:
list_anti_u2 = []
list_anti_ne2 = []
with open('ANTI_U2.txt','r') as f:
    fu = [line[:-1] for line in f]

for i in range(len(fu)):
    list_anti_u2.append(float(fu[i]))

with open('ANTI_ne2.txt','r') as f:
    fne = [line[:-1] for line in f]

for i in range(len(fu)):
    list_anti_ne2.append(float(fne[i]))

anti_tu2 = []
anti_ne2 = []
for j in range(len(fu)):
    if (list_anti_u2[j]!=-1 and list_anti_ne2[j]!=-1):
        anti_tu2.append(t/list_anti_u2[j])
        anti_ne2.append(list_anti_ne2[j])
########
#       1º region:
list_para_u1 = []
list_para_ne1 = []
with open('PARA_U1.txt','r') as f:
    fu = [line[:-1] for line in f]

for i in range(len(fu)):
    list_para_u1.append(float(fu[i]))

with open('PARA_ne1.txt','r') as f:
    fne = [line[:-1] for line in f]

for i in range(len(fu)):
    list_para_ne1.append(float(fne[i]))

para_tu1 = []
para_ne1 = []
for j in range(len(fu)):
    if (list_para_u1[j]!=-1 and list_para_ne1[j]!=-1):
        para_tu1.append(t/list_para_u1[j])
        para_ne1.append(list_para_ne1[j])

#       2º region:
list_para_u2 = []
list_para_ne2 = []
with open('PARA_U2.txt','r') as f:
    fu = [line[:-1] for line in f]

for i in range(len(fu)):
    list_para_u2.append(float(fu[i]))

with open('PARA_ne2.txt','r') as f:
    fne = [line[:-1] for line in f]

for i in range(len(fu)):
    list_para_ne2.append(float(fne[i]))

para_tu2 = []
para_ne2 = []
for j in range(len(fu)):
    if (list_para_u2[j]!=-1 and list_para_ne2[j]!=-1):
        para_tu2.append(t/list_para_u2[j])
        para_ne2.append(list_para_ne2[j])
########
label = [0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8]
str_label=list(map(str,label))
string = '$'
str_label = [string + x + string for x in str_label]

ax = plt.subplot()
#ax.set_ylim(ymin=None, ymax=0.215)
ax.scatter(anti_ne1, anti_tu1, color="blue")
ax.scatter(para_ne1, para_tu1, color="lime")
ax.scatter(ferro_ne1, ferro_tu1, color="red")
ax.scatter(ferro_ne2, ferro_tu2, color="red")
ax.scatter(anti_ne2, anti_tu2, color="blue")
ax.scatter(para_ne2, para_tu2, color="lime")
ax.set_ylabel(r'$t/U$')
ax.set_xlabel(r'$n_e/2$')
ax.set_xticks(label)
ax.set_xticklabels(str_label)
ant = plt.text(0.39, 0.76, 'Antiferromagnetic', transform=ax.transAxes, fontsize=12)
ant.set_bbox(dict(facecolor='white', alpha=0.93, edgecolor='white'))
fer1 = plt.text(0.17, 0.32, 'Ferromagnetic', transform=ax.transAxes, fontsize=12)
fer1.set_bbox(dict(facecolor='white', alpha=0.95, edgecolor='white'))
fer2 = plt.text(0.7, 0.32, 'Ferromagnetic', transform=ax.transAxes, fontsize=12)
fer2.set_bbox(dict(facecolor='white', alpha=0.95, edgecolor='white'))
par1 = plt.text(0.07, 0.86, 'Paramagnetic', transform=ax.transAxes, fontsize=12)
par1.set_bbox(dict(facecolor='white', alpha=0.9, edgecolor='white'))
par2 = plt.text(0.75, 0.86, 'Paramagnetic', transform=ax.transAxes, fontsize=12)
par2.set_bbox(dict(facecolor='white', alpha=0.9, edgecolor='white'))
#ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),
#          ncol=3, fancybox=True, shadow=True)

plt.show()
