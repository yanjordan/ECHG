import sys
import os
import string
import numpy as np 
import numpy.matlib
from matplotlib import pyplot as plt
from PyQt5 import uic, QtCore, QtGui, QtWidgets
from ECHG import Ui_MainWindow
from PyQt5.QtWidgets import QFileDialog

Ui_MainWindow, QtBaseClass = uic.loadUiType("ECHG.ui")

""" class MainWindow(QtWidgets.QWidget):
    def __init__(self):
        super(MainWindow,self).__init__()

        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        self.ui.Select_exe.clicked.connect(self.exe_choose) """
class MyECHG(QtWidgets.QMainWindow, Ui_MainWindow):

    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)
        self.exe_path.setText(os.getcwd() + os.sep + 'bin' + os.sep + 'get_ECHG.exe')

        self.Select_exe.clicked.connect(self.exe_choose)
        self.Select_dossier.clicked.connect(self.dossier_choose)
        self.Select_ECHG.clicked.connect(self.ECHG_choose)
        self.Select_PATO.clicked.connect(self.PATO_choose)
        self.Select_Atom.clicked.connect(self.Atom_choose)
        self.Run.clicked.connect(self.Run_function)



    def exe_choose(self):
        filename,_ = QFileDialog.getOpenFileName(self,'exe file','.',"exe file (*.exe)")
        self.exe_path.setText(filename)
        
    def dossier_choose(self):
        dossier_path = QFileDialog.getExistingDirectory(self,'data folder', '.')
        self.Dossier_path.setText(dossier_path)

    def ECHG_choose(self):
        filename,_ = QFileDialog.getOpenFileName(self,'ECHG data file',self.Dossier_path.toPlainText(),"exe file (*.data *.dat)")
        self.ECHG_path.setText(filename)

    def PATO_choose(self):
        filename,_ = QFileDialog.getOpenFileName(self,'PATO data file',self.Dossier_path.toPlainText(),"exe file (*.data *.dat)")
        self.PATO_path.setText(filename)

    def Atom_choose(self):
        filename,_ = QFileDialog.getOpenFileName(self,'Atom list file',self.Dossier_path.toPlainText(),"exe file (*.txt)")
        self.Atom_path.setText(filename)

    def Run_function(self):
        if self.user_define.isChecked():
            intel_val = float(self.User_val.toPlainText())
            level_pos = np.arange(intel_val, intel_val*20, intel_val)
            level_neg = -level_pos
            level_neg.sort()
        else:
            intel_val = float(self.Bader_val.toPlainText())
            level_pos = intel_val*2**np.arange(0, 20, 1)
            level_neg = -level_pos
            level_neg.sort()

        Lwidth = float(self.Line_width.toPlainText())
        Hlength = float(self.H_bond.toPlainText())
        Blength = float(self.Bond.toPlainText())

        Cwidth = float(self.contour_width.toPlainText())
        Fsize = int(self.fontsize.toPlainText())
        Asize = int(self.axiesize.toPlainText())
        Reso = int(self.resolution.toPlainText())
        Pdis = float(self.plan_dis.toPlainText())

        Pcolor = self.pos_color.currentText()
        Pline = self.pos_line.currentText().split(" ")[0]

        Ncolor = self.neg_color.currentText()
        Nline = self.neg_line.currentText().split(" ")[0]

        Necolor = self.neu_color.currentText()
        Neline = self.neu_line.currentText().split(" ")[0]

        #print([Pline, Nline, Neline])

        Element=['H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe']

        os.system(self.exe_path.toPlainText() + ' ' + self.ECHG_path.toPlainText())

        ##########  data plot ######
        unit = 0.529177**3
        
        Atome = np.loadtxt('Atome.dat')
        Lattice = np.loadtxt('Lattice.dat')
        Parameters = np.loadtxt('Parameters.dat')

        #### atomes list
        

        AtomeDict = {}
        if os.path.exists(self.Atom_path.toPlainText()):
            dicFile = open(self.Atom_path.toPlainText(), 'r')

            while True:
                line = dicFile.readline()
                if line == '':
                    break
                index = line.find('\t')
                key = line[:index]
                value = line[index:]

                key = key.replace('\t','').replace('\n','')
                key_number = int(key)
                value = value.replace('\t','').replace('\n','')
                AtomeDict[key_number] = value

            dicFile.close()
        else:
            print('There is no atome list file.')

        Ihferm = int(Parameters[0])
        Nrow = int(Parameters[1])
        Ncol = int(Parameters[2])
        Dx = Parameters[3]
        Dy = Parameters[4]
        Cosxy = Parameters[5]
        A = Parameters[6:9]
        B = Parameters[9:12]
        C = Parameters[12:15]
        Naf = int(Parameters[15])
        Ldim = int(Parameters[16])
        #print([Ihferm, Nrow, Ncol, Naf, Ldim])
        Atomes = np.zeros((Naf*27, 4))

        Density = np.loadtxt('Density.dat')
        if (Ihferm ==1 or Ihferm ==3):
            Spin = np.loadtxt('Spin.dat')
        else:
            Spin = np.zeros((Nrow, Ncol))

        ### corrected the atomic number which use the ECP basis set
        for i in range(Naf):
            if Atome[i, 0] > 300:
                Atome[i, 0] = Atome[i, 0] - 300
            elif Atome[i, 0] > 200:
                Atome[i, 0] = Atome[i, 0] - 200

        #### generate larger atomic list of le supercell (27 cells)
        n_mat=0
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    deplace = (i-1)*Lattice[0:1,:] + (j-1)*Lattice[1:2,:] + (k-1)*Lattice[2:3,:]
                    Atomes[n_mat*Naf:(n_mat+1)*Naf,0:1]=Atome[0:Naf,0:1]
                    Atomes[n_mat*Naf:(n_mat+1)*Naf,1:4]=Atome[0:Naf,1:4]+np.matlib.repmat(deplace, Naf, 1)
                    n_mat=n_mat+1
        
        os.system(self.exe_path.toPlainText() + ' ' + self.PATO_path.toPlainText())
        Den_PATO = np.loadtxt('Density.dat')

        Density = Density/unit
        Spin = Spin/unit
        Den_PATO = Den_PATO/unit

        Defor = Density - Den_PATO

        #np.savetxt('Atom.out', Atome, fmt='%4d % 14.12e % 14.12e % 14.12e')
        #np.savetxt('Atomes.out', Atomes, fmt='%4d % 14.12e % 14.12e % 14.12e')
        
        if (Ihferm ==1 or Ihferm ==3):
            os.remove('Spin.dat')

        os.remove('Atome.dat')
        os.remove('Lattice.dat')
        os.remove('Parameters.dat')
        os.remove('Density.dat')
        
        #####
        if (Ihferm ==1 or Ihferm ==3):
            Up = (Density + Spin)/2 - Den_PATO/2
            Down = (Density - Spin)/2 - Den_PATO/2
        else:
            Up = Density/2 - Den_PATO/2
            Down = Density/2 - Den_PATO/2

        D_max = np.max(Defor)
        D_min = np.min(Defor)

        S_max = np.max(Spin)
        S_min = np.min(Spin)

        Do_max = np.max(Down)
        Do_min = np.min(Down)

        U_max = np.max(Up)
        U_min = np.min(Up)

        ### define vectors of plane

        Vector_x = A - B
        Vector_x = Vector_x/np.linalg.norm(Vector_x)

        Vector_y = C - B
        Vector_y = Vector_y/np.linalg.norm(Vector_y)

        Vector_z = np.cross(Vector_x, Vector_y)

        ### project the corrdinates into the plane

        for i in range(Naf*27):
            Vector = Atomes[i, 1:4] - B 
            Atomes[i,1:2] = np.dot(Vector, Vector_y)
            Atomes[i,2:3] = np.dot(Vector, Vector_x)
            Atomes[i,3:4] = np.dot(Vector, Vector_z)

        X_max = np.linalg.norm(A-B)
        Y_max = np.linalg.norm(C-B)

        ### define contour inter-value
        if self.user_define.isChecked():
            interval = float(self.User_val.toPlainText())
            D_pos = np.arange(interval, interval*20, interval)
            D_neg = -D_pos
            D_neg.sort()

            S_pos = np.arange(interval, interval*20, interval)
            S_neg = -S_pos
            S_neg.sort()

            Do_pos = np.arange(interval, interval*20, interval)
            Do_neg = -Do_pos
            Do_neg.sort()

            U_pos = np.arange(interval, interval*20, interval)
            U_neg = -U_pos
            U_neg.sort()

        elif self.Bader.isChecked():
            interval = float(self.Bader_val.toPlainText())
            D_pos = interval*2**np.arange(0, 20, 1)
            D_neg = -D_pos
            D_neg.sort()

            S_pos = interval*2**np.arange(0, 20, 1)
            S_neg = -S_pos
            S_neg.sort()

            Do_pos = interval*2**np.arange(0, 20, 1)
            Do_neg = -Do_pos
            Do_neg.sort()

            U_pos = interval*2**np.arange(0, 20, 1)
            U_neg = -U_pos
            U_neg.sort()

        ### define plot parameters
        x = np.linspace(0, Dx*(Nrow-1), Nrow)
        y = np.linspace(0, Dy*(Ncol-1), Ncol)
        X, Y = np.meshgrid(y, x)

        if self.Defo.isChecked():
            plt.contour(X,Y,Defor[::-1], levels = D_pos, colors = Pcolor, linestyles = Pline, linewidths = Cwidth)
            plt.contour(X,Y,Defor[::-1], levels = D_neg, colors = Ncolor, linestyles = Nline, linewidths = Cwidth)
            plt.contour(X,Y,Defor[::-1], levels = [0.0], colors = Necolor, linestyles = Neline, linewidths = Cwidth)
        elif self.Spin.isChecked():
            plt.contour(X,Y,Spin[::-1], levels = S_pos, colors = Pcolor, linestyles = Pline, linewidths = Cwidth)
            plt.contour(X,Y,Spin[::-1], levels = S_neg, colors = Ncolor, linestyles = Nline, linewidths = Cwidth)
            plt.contour(X,Y,Spin[::-1], levels = [0.0], colors = Necolor, linestyles = Neline, linewidths = Cwidth)
        elif self.Up.isChecked():
            plt.contour(X,Y,Up[::-1], levels = U_pos, colors = Pcolor, linestyles = Pline, linewidths = Cwidth)
            plt.contour(X,Y,Up[::-1], levels = U_neg, colors = Ncolor, linestyles = Nline, linewidths = Cwidth)
            plt.contour(X,Y,Up[::-1], levels = [0.0], colors = Necolor, linestyles = Neline, linewidths = Cwidth)
        elif self.Down.isChecked():
            plt.contour(X,Y,Down[::-1], levels = Do_pos, colors = Pcolor, linestyles = Pline, linewidths = Cwidth)
            plt.contour(X,Y,Down[::-1], levels = Do_neg, colors = Ncolor, linestyles = Nline, linewidths = Cwidth)
            plt.contour(X,Y,Down[::-1], levels = [0.0], colors = Necolor, linestyles = Neline, linewidths = Cwidth)

        if not self.notitle.isChecked():
            if self.user_define.isChecked():
                plt.title('Iso-contour are every ' + self.User_val.toPlainText() + 'e/A^3', fontsize = Fsize)
            elif self.Bader.isChecked():
                plt.title('Iso-contour of Bader: $\pm$' + self.Bader_val.toPlainText() + 'X 2^n e/A^3', fontsize = Fsize)
        ### draw atomes
        Atomes_tag = np.zeros(Naf*27, dtype=int)

        if not self.noatoms.isChecked():
            N_tag=0
            for i in range(Naf*27):
                if (Atomes[i, 1:2] >= 0.0 and Atomes[i, 1:2] <= Y_max) :
                    if (Atomes[i, 2:3] >= 0.0 and Atomes[i, 2:3] <= X_max) :
                        if abs(Atomes[i, 3:4]) < Pdis :
                            Atomes_tag[N_tag] = i
                            N_tag = N_tag+1
                            
                            atome_num = i%Naf+1
                            atomic_num = int(Atome[atome_num,0:1])
       

                            if abs(Atomes[i,3:4]) <0.001 :
                                plt.text(Atomes[i,1:2], Atomes[i,2:3], Element[atomic_num-1] + AtomeDict.get(atome_num, str(atome_num)), fontsize = Fsize)
                                plt.scatter(Atomes[i,1:2], Atomes[i,2:3], marker = '.')
                            elif Atomes[i,3:4] < 0.0 :
                                plt.text(Atomes[i,1:2], Atomes[i,2:3], Element[atomic_num-1] + AtomeDict.get(atome_num, str(atome_num)), fontsize = Fsize)
                                plt.scatter(Atomes[i,1:2], Atomes[i,2:3], marker = '*')
                            else:
                                plt.text(Atomes[i,1:2], Atomes[i,2:3], Element[atomic_num-1] + AtomeDict.get(atome_num, str(atome_num)), fontsize = Fsize)
                                plt.scatter(Atomes[i,1:2], Atomes[i,2:3], marker = 'x')
                            
            ### draw bonds
            line_x = np.zeros(2)
            line_y = np.zeros(2)
            for i in range(N_tag):
                for j in range(i+1,N_tag+1):
                    dis_atomes = np.linalg.norm(Atomes[Atomes_tag[i]:Atomes_tag[i]+1, 1:3] - Atomes[Atomes_tag[j]:Atomes_tag[j]+1, 1:3])

                    line_x[0] = Atomes[Atomes_tag[i]:Atomes_tag[i]+1 ,1:2]
                    line_x[1] = Atomes[Atomes_tag[j]:Atomes_tag[j]+1 ,1:2]
                    line_y[0] = Atomes[Atomes_tag[i]:Atomes_tag[i]+1 ,2:3]
                    line_y[1] = Atomes[Atomes_tag[j]:Atomes_tag[j]+1 ,2:3]

                    if (Atomes[Atomes_tag[i]:Atomes_tag[i]+1, 0:1] == 1 or Atomes[Atomes_tag[j]:Atomes_tag[j]+1, 0:1] == 1):
                        if dis_atomes < Hlength :
                            plt.plot(line_x, line_y,'-', color = 'Black', linewidth = Lwidth)
                    else:
                        if dis_atomes < Blength :
                            plt.plot(line_x, line_y,'-', color = 'Black', linewidth = Lwidth)
        plt.axis('equal')
        #plt.xlim(0.0, X_max)
        #plt.ylim(0.0, Y_max)

        ### correct the figure
        if self.up2down.isChecked():
            plt.gca().invert_yaxis()
        if self.left2right.isChecked():
            plt.gca().invert_xaxis()

        if not self.xyticks.isChecked():
            plt.xticks([])
            plt.yticks([])

        if self.boxwhite.isChecked():
            plt.axis('off')

        

        if self.Defo.isChecked():
            if self.ps.isChecked():
                plt.savefig(self.Dossier_path.toPlainText() + '/Defromation.ps', dpi =  Reso)
            if self.png.isChecked():
                plt.savefig(self.Dossier_path.toPlainText() + '/Defromation.png', dpi =  Reso)
            if self.eps.isChecked():
                plt.savefig(self.Dossier_path.toPlainText() + '/Defromation.eps', dpi =  Reso)
            if self.pdf.isChecked():
                plt.savefig(self.Dossier_path.toPlainText() + '/Defromation.pdf', dpi =  Reso)
        elif self.Spin.isChecked():
            if self.ps.isChecked():
                plt.savefig(self.Dossier_path.toPlainText() + '/Spin.ps', dpi =  Reso)
            if self.png.isChecked():
                plt.savefig(self.Dossier_path.toPlainText() + '/Spin.png', dpi =  Reso)
            if self.eps.isChecked():
                plt.savefig(self.Dossier_path.toPlainText() + '/Spin.eps', dpi =  Reso)
            if self.pdf.isChecked():
                plt.savefig(self.Dossier_path.toPlainText() + '/Spin.pdf', dpi =  Reso)    
        elif self.Up.isChecked():
            if self.ps.isChecked():
                plt.savefig(self.Dossier_path.toPlainText() + '/Up.ps', dpi =  Reso)
            if self.png.isChecked():
                plt.savefig(self.Dossier_path.toPlainText() + '/Up.png', dpi =  Reso)
            if self.eps.isChecked():
                plt.savefig(self.Dossier_path.toPlainText() + '/Up.eps', dpi =  Reso)
            if self.pdf.isChecked():
                plt.savefig(self.Dossier_path.toPlainText() + '/Up.pdf', dpi =  Reso)
        elif self.Down.isChecked():
            if self.ps.isChecked():
                plt.savefig(self.Dossier_path.toPlainText() + '/Down.ps', dpi =  Reso)
            if self.png.isChecked():
                plt.savefig(self.Dossier_path.toPlainText() + '/Down.png', dpi =  Reso)
            if self.eps.isChecked():
                plt.savefig(self.Dossier_path.toPlainText() + '/Down.eps', dpi =  Reso)
            if self.pdf.isChecked():
                plt.savefig(self.Dossier_path.toPlainText() + '/Down.pdf', dpi =  Reso)
        
        
        plt.show()


if __name__ == "__main__":
    app=QtWidgets.QApplication.instance()
    if not app:
        app = QtWidgets.QApplication(sys.argv)

    window = MyECHG()
    window.show()
    sys.exit(app.exec_())

""" app = QtWidgets.QApplication(sys.argv)
MainWindow = QtWidgets.QMainWindow()
ui = Ui_MainWindow()
ui.setupUi(MainWindow)
MainWindow.show()
sys.exit(app.exec_()) """
