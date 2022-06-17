import matplotlib.pyplot as plt
import csv

x = []
y_gamma = []
y_mu_minus = []
y_mu_plus = []
y_nu_e = []
y_nu_mu = []
y_anti_nu_e = []
y_anti_nu_mu = []

#choose the graphs you want putting true:
gamma_graph= True
mu_graph = True
neutrino_graph = True

with open("gamma.txt", "r") as csvfile:
    plots = csv.reader(csvfile, delimiter=" ")
    for row in plots:
        y_gamma.append(float(row[0]))
        x.append(float(row[1]))

with open("mu+.txt", "r") as csvfile:
    plots = csv.reader(csvfile, delimiter=" ")
    for row in plots:
        y_mu_plus.append(float(row[0]))


with open("mu-.txt", "r") as csvfile:
    plots = csv.reader(csvfile, delimiter=" ")
    for row in plots:
        y_mu_minus.append(float(row[0]))

with open("nu_e.txt", "r") as csvfile:
    plots = csv.reader(csvfile, delimiter=" ")
    for row in plots:
        y_nu_e.append(float(row[0]))

with open("nu_mu.txt", "r") as csvfile:
    plots = csv.reader(csvfile, delimiter=" ")
    for row in plots:
        y_nu_mu.append(float(row[0]))

with open("anti_nu_e.txt", "r") as csvfile:
    plots = csv.reader(csvfile, delimiter=" ")
    for row in plots:
        y_anti_nu_e.append(float(row[0]))

with open("anti_nu_mu.txt", "r") as csvfile:
    plots = csv.reader(csvfile, delimiter=" ")
    for row in plots:
        y_anti_nu_mu.append(float(row[0]))



if (gamma_graph):
    #gamma
    plt.scatter(x, y_gamma, s=20)
    plt.title("Geant4 hadronic shower graphs")
    plt.ylabel("Nº gamma")
    plt.xlabel("distance (km)")
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    plt.show()

if (mu_graph):
    #muons 
    line1 = plt.scatter(x, y_mu_plus, c='b', label='1', s=20, marker = '^')
    line2 = plt.scatter(x,y_mu_minus, c='r', label='1', s=20)
    plt.title("Geant4 hadronic shower graphs")
    plt.ylabel("Nº muons")
    plt.xlabel("distance (km)")
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    plt.legend((line1, line2), ('mu-', 'mu+'))
    plt.show()

if (neutrino_graph):
    #neutrinos
    line1 = plt.scatter(x, y_nu_e, c='b', label='1', s=20)
    line2 = plt.scatter(x,y_anti_nu_e, c='r', label='1', s=20)
    plt.title("Geant4 hadronic shower graphs")
    plt.ylabel("Nº electron neutrinos")
    plt.xlabel("distance (km)")
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    plt.legend((line1, line2), ('nu_e', 'anti_nu_e'))
    plt.show()

    line1 = plt.scatter(x, y_nu_mu, c='b', label='1', s=20)
    line2 = plt.scatter(x,y_anti_nu_mu, c='r', label='1', s=20, marker = '^')
    plt.title("Geant4 hadronic shower graphs")
    plt.ylabel("Nº muon neutrinos")
    plt.xlabel("distance (km)")
    plt.xlim(xmin=0)
    plt.ylim(ymin=0)
    plt.legend((line1, line2), ('nu_mu', 'anti_nu_mu'))
    plt.show()
