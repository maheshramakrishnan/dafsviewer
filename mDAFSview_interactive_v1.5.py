# -*- coding: utf-8 -*-
"""
Development of interactive version of the integration SW started Nov 2022
@author: mahesh.ramakrishnan@hotmail.com

Shows live data and integrated images side by side.

Version history: 
    + v1.0: Developmental tool  
    + v1.1: Trying to add more panels to plot other things
    + v1.2: Adding markers for slider update
    + v1.3: Many features added.
    + v1.4: More features like peak selection.
    + v1.5: I0, I1 in separate panel, better arrangement of panels in layout.
            
To do (wishlist):
    - Panels for I0, I1, mu, DAFS
    - Position markers (Energy) with slider update
    - Energy value on slider
    - Mouse hover point readout for all image windows.
    - Fitting info for the current curve as a label widget (not in console).
    - Peakshape asymmetry.
"""
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets
import pyqtgraph.console
import numpy as np
import time, os, h5py, hdf5plugin
from scipy.optimize import curve_fit
from scipy.signal import savgol_filter

# List of global variables
result_file = ''
alba_file = ''
int_data = []
peak1_I = []
energy = []
I0 = []
I1 = []
Int_time = []
mu = []
N_saved_curves = np.zeros((1),dtype='uint32')
q = np.zeros((3000))

# Opening QT application and defining GUI    
app = QtGui.QApplication([])
pg.setConfigOption('leftButtonPan', False)

# Set up UI widgets
win = pg.QtGui.QWidget()
win.setWindowTitle('DAFSViewer_v1.3')
icon = QtGui.QIcon('dafs.png')
win.setWindowIcon(icon)
layout = pg.QtGui.QGridLayout()

for i in range(3,9):
    layout.setRowMinimumHeight(i,80)
   
win.setLayout(layout)
layout.setContentsMargins(0, 0, 0, 0)
pg.setConfigOption('imageAxisOrder', 'row-major')


# Defining slider widget, image number display label

sliderLabel = pg.QtGui.QLabel('Energy selector:')
sliderLabel.setAlignment(QtCore.Qt.AlignRight)
layout.addWidget(sliderLabel, 2, 0)
slider = pg.Qt.QtGui.QSlider(QtCore.Qt.Horizontal)
slider.setValue(0)
layout.addWidget(slider, 2, 1, 1, 4)

sliderVal = pg.QtGui.QLabel(str(slider.value()))
sliderValFont = sliderVal.font()
sliderValFont.setBold(True)
sliderValFont.setPointSize(10)
sliderVal.setFont(sliderValFont)
layout.addWidget(sliderVal, 2, 5, 1, 2)

fitVal = pg.QtGui.QLabel('PeakFit parameters appear here.')
layout.addWidget(fitVal, 6, 2, 2, 1)

minLabel = pg.QtGui.QLabel('Peak Center min:')
minLabel.setAlignment(QtCore.Qt.AlignRight)
layout.addWidget(minLabel, 9, 1)
minSpin = pg.SpinBox(value=1845, step=1, bounds=[0, 2999], delay=0, int=True)
layout.addWidget(minSpin, 9, 2)

maxLabel = pg.QtGui.QLabel('Peak Center max:')
maxLabel.setAlignment(QtCore.Qt.AlignRight)
layout.addWidget(maxLabel, 10, 1)
maxSpin = pg.SpinBox(value=2185, step=1, bounds=[0, 2999], delay=0, int=True)
layout.addWidget(maxSpin, 10, 2)

a1 = pg.ArrowItem(angle=180, tipAngle=60, headLen=5, tailLen=20, tailWidth=5, pen={'color': 'k', 'width': 3})
a3 = pg.ArrowItem(angle=90, tipAngle=30, headLen=5, tailLen=20, tailWidth=5, pen={'color': 'k', 'width': 2})
a4 = pg.ArrowItem(angle=90, tipAngle=30, headLen=5, tailLen=20, tailWidth=5, pen={'color': 'k', 'width': 2})
a51 = pg.ArrowItem(angle=90, tipAngle=30, headLen=5, tailLen=20, tailWidth=5, pen={'color': 'r', 'width': 2})
a52 = pg.ArrowItem(angle=90, tipAngle=30, headLen=5, tailLen=20, tailWidth=5, pen={'color': 'b', 'width': 2})


# Defining text widgets for path, filenames
resultLabel = pg.QtGui.QLabel('Result File:')
resultLabel.setAlignment(QtCore.Qt.AlignRight)
layout.addWidget(resultLabel, 0, 0)
resultLine = pg.Qt.QtGui.QLineEdit('C:\\Users\\mahram\\Desktop\\data\DAFS\\perovskite\\202212\\MAPbBrI_2ME_DMSO_coat015a_DAFS.h5')
layout.addWidget(resultLine, 0, 1, 1, 2)

albaLabel = pg.QtGui.QLabel('Alba (sardana master) File:')
albaLabel.setAlignment(QtCore.Qt.AlignRight)
layout.addWidget(albaLabel, 1, 0)
albaLine = pg.Qt.QtGui.QLineEdit('C:\\Users\\mahram\\Desktop\\data\DAFS\\perovskite\\202212\\MAPbBrI_2ME_DMSO_coat015a_DAFS_resultCluster.h5')
layout.addWidget(albaLine, 1, 1, 1, 2)

# Buttons and check-boxes

loadResButton = QtGui.QPushButton('Load Result File')
loadResButton.resize(80,20)
layout.addWidget(loadResButton,0,3)

calcButton = QtGui.QPushButton('Generate DAFS')
calcButton.resize(80,20)
layout.addWidget(calcButton,0,4)

progBar = QtWidgets.QProgressBar()
progBar.setMaximumHeight(20)
layout.addWidget(progBar,0,5)

saveCurveButton = QtGui.QPushButton('Save Current Curve')
layout.addWidget(saveCurveButton,9,0)

clrCurveButton = QtGui.QPushButton('Clear All Curves')
layout.addWidget(clrCurveButton,10,0)

logCheckBox = QtGui.QCheckBox('Integrated Data Log Scale')
# logCheckBox.clicked(1)
layout.addWidget(logCheckBox,1,3)

normCheckBox = QtGui.QCheckBox('Normalize I0')
layout.addWidget(normCheckBox,1,4)

normCheckBox2 = QtGui.QCheckBox('Normalize dt')
layout.addWidget(normCheckBox2,1,5)

absCheckBox = QtGui.QCheckBox('Abs. Corr. DAFS')
layout.addWidget(absCheckBox,12,1)

savgolCheckBox = QtGui.QCheckBox('Savitzky-Golay Filter')
layout.addWidget(savgolCheckBox,12,2)

# Setting the console widget
consInitText = 'Välkommen till sveriges bästa DAFSViewer!\n'
cons = pyqtgraph.console.ConsoleWidget(text=consInitText)
layout.addWidget(cons, 15, 0, 1, 3)

# Savgol filter selection

savgolWinLabel = pg.QtGui.QLabel('Sav-Gol Window length:')
savgolWinLabel.setAlignment(QtCore.Qt.AlignRight)
layout.addWidget(savgolWinLabel, 13, 1)
savgolWin = pg.ComboBox()
savgolWin.setItems({'3':3,'7':7,'11':11,'15':15})
savgolWin.setValue(15)
savgolWin.resize(40, 20)
layout.addWidget(savgolWin, 13, 2)

savgolPolLabel = pg.QtGui.QLabel('Sav-Gol Polynom order:')
savgolPolLabel.setAlignment(QtCore.Qt.AlignRight)
layout.addWidget(savgolPolLabel, 14, 1)
savgolPol = pg.ComboBox()
savgolPol.setItems({'1':1,'2':2,'3':3})
savgolPol.setValue(3)
savgolPol.resize(40, 20)
layout.addWidget(savgolPol, 14, 2)

# Setting the plot widgets

w1 = pg.GraphicsLayoutWidget()
layout.addWidget(w1, 3, 0, 3, 3)
w1.resize(1500, 1500)
# w1.ci.setBorder((50, 50, 100))
w1.setBackground((255,250,250))
p1 = w1.addPlot(title="Azimuthal Integrated Data")
img1 = pg.ImageItem()
p1.addItem(img1)
hist1 = pg.HistogramLUTItem()
hist1.setImageItem(img1)
w1.addItem(hist1)
p1.setLabel('bottom',text="Radial Bins")
p1.setLabel('left',text="Energy point number")

w2 = pg.GraphicsLayoutWidget()
layout.addWidget(w2, 6, 0, 3, 2)
w2.resize(1000, 1500)
# w2.ci.setBorder((50, 50, 100))
w2.setBackground((255,250,250))
p2 = w2.addPlot(title="Line Plots")
p2.setLabel('bottom',text="Radial axis")
#p2.setAspectLocked(1.0)

w3 = pg.GraphicsLayoutWidget()
w3.setBackground((255,250,250))
w3.resize(1000,1500)
layout.addWidget(w3, 6, 3, 3, 3)
p3 = w3.addPlot(title="mu")
p3.setLabel('bottom',text="Energy (eV)")
#p3.setAspectLocked(1.0)

w4 = pg.GraphicsLayoutWidget()
w4.setBackground((255,250,250))
w4.resize(1000,1500)
layout.addWidget(w4, 9, 3, 7, 3)
p4 = w4.addPlot(title="DAFS")
p4.setLabel('bottom',text="Energy (eV)")

w5 = pg.GraphicsLayoutWidget()
w5.setBackground((255,250,250))
w5.resize(1000,1500)
layout.addWidget(w5, 3, 3, 3, 3)
p5 = w5.addPlot(title="Raw currents")
p5.setLabel('bottom',text="Energy (eV)")



# Peak functions
def f_gauss(x,A,B,C,D,E):
    y = A*np.exp(np.divide(-1*np.square(x-B),2*np.square(C)))+D*x+E
#    y = A*np.exp(np.divide(-1*np.square(x-B),2*np.square(C)))+E
    return y

def f_psVgt(x,A1,A2,B,C1,C2,D,E):
    y = A1*np.exp(np.divide(-np.square(x-B),2*np.square(C1))) + A2*np.divide(C2,(np.square(x-B)+np.square(C2)))+D*x+E
    return y
   

# User input and data loading

def load_textlines():
    global result_file, alba_file
    if resultLine.text():
        result_file = resultLine.text()
    if albaLine.text():
        alba_file = albaLine.text()

def load_resultData():
    global int_data, q, img1, I0, I1, Int_time, mu, energy
    load_textlines()
    cons.write('Loading result image, please wait..\n')
    
    with h5py.File(result_file, 'r') as fh:
        rawresult_full = np.array(fh['/Rawresult'])
        q = np.array(fh['/q (E = 7601.07 eV)'])
    
    with h5py.File(alba_file, 'r') as fh:
        I0 = np.array(fh['/I0a'])+np.array(fh['/I0b'])
        I1 = np.array(fh['/I1'])+np.array(fh['/It'])
        Int_time = np.array(fh['/integration_time'])
        energy = np.array(fh['/energy'])
    
    energy = energy[0]
    mu = np.log(I0[0]/I1[0])
    
    raw0 = np.split(rawresult_full,6)[0]
    raw1 = np.split(rawresult_full,6)[1]
    raw2 = np.split(rawresult_full,6)[2]
    raw3 = np.split(rawresult_full,6)[3]
    raw4 = np.split(rawresult_full,6)[4]
    raw5 = np.split(rawresult_full,6)[5]
    int_data = raw0+raw1+raw2+raw3+raw4+raw5
    
    if normCheckBox.checkState():
        cons.write('Normalizing with I0..')
        for i in range(0,len(I0[0])):
            int_data[i] = (raw0[i]/I0[0,i]+raw1[i]/I0[1,i]+raw2[i]/I0[2,i]+raw3[i]/I0[3,i]+raw4[i]/I0[4,i]+raw5[i]/I0[5,i])*0.000001
    
    if normCheckBox2.checkState():
        cons.write('Normalizing with It (intg. time)..')
        for i in range(0,len(Int_time[0])):
            int_data[i] = (raw0[i]/Int_time[0,i]+raw1[i]/Int_time[1,i]+raw2[i]/Int_time[2,i]+raw3[i]/Int_time[3,i]+raw4[i]/Int_time[4,i]+raw5[i]/Int_time[5,i])*1000000
    
    log_resultData()

  
def log_resultData():        
    global int_data, img1, slider, p3, mu, I0, I1
    if logCheckBox.checkState():
        img1.setImage(np.log10(int_data))
        hist1.setLevels(0.01*np.max(np.log10(int_data)), 0.8*np.max(np.log10(int_data)))
    else:
        img1.setImage(int_data)
        hist1.setLevels(0.1*np.max(int_data), 0.8*np.max(int_data))
    img1.setOpts(colorMap='inferno')
    p1.show()
    p1.setXRange(0,2999)
    
    p3.plot(energy, mu, pen=(255,0,0))
    p3.setXRange(energy[0],energy[int(len(energy)/2)])
    p3.setYRange(mu.min(),mu.max())
    
    p3.show()
    p5.plot(energy, np.mean(I0, axis=0), pen=(255,0,0))
    p5.plot(energy, np.mean(I1,axis=0), pen=(0,0,255))
    p5.show()
    
    slider.setMinimum(0)
    slider.setMaximum(len(int_data)-1)
    slider_update()
    cons.write('Results Loaded.\n')

def fit_PseudoVoigt_all():
    global int_data, peak1_I, p3, p4, energy, mu
    cons.write('Fitting DAFS. Status displayed every 100 energy points..')
    peak1_cen = int((minSpin.value()+maxSpin.value())/2)
    peak1_disp_up = np.linspace(maxSpin.value(),minSpin.value(),3599, dtype='float32')
    peak1_disp_down = np.flip(peak1_disp_up)
    peak1_disp = np.append(peak1_disp_up,peak1_disp_down)
    
    for i in range(0,len(int_data)):
        x1 = np.linspace(int(peak1_disp[i])-200,int(peak1_disp[i])+199,400)
        y1 = int_data[i][int(peak1_disp[i])-200:int(peak1_disp[i])+200]
        popt, pcov = curve_fit(f_psVgt, x1, y1, bounds=([0.1,0.1,peak1_disp[i]-50,0.5,0.5,-100,-1000],[1000000,1000000,peak1_disp[i]+50,30,30,100,5000]))
    
        I = f_psVgt(x1,popt[0],popt[1],popt[2],popt[3],popt[4],0,0)
        peak1_I.append(np.trapz(I,x1))
        progBar.setValue((i/(len(int_data)-1))*100)
        #perr = np.sqrt(np.diag(pcov))
        #cons.write('Standard deviation:')
        #cons.write(perr)  
    p3.plot(energy, mu, pen=(255,0,0))
    p4.plot(energy, peak1_I, pen=(0,0,255))
    p3.show()
    p4.show()
        
def fit_PseudoVoigt_single(i):
    global int_data, p2
    
    peak1_cen = int((minSpin.value()+maxSpin.value())/2)
    peak1_disp_up = np.linspace(maxSpin.value(),minSpin.value(),3599, dtype='float32')
    peak1_disp_down = np.flip(peak1_disp_up)
    peak1_disp = np.append(peak1_disp_up,peak1_disp_down)
    
    x1 = np.linspace(int(peak1_disp[i])-200,int(peak1_disp[i])+199,400)
    y1 = int_data[i][int(peak1_disp[i])-200:int(peak1_disp[i])+200]
    popt, pcov = curve_fit(f_psVgt, x1, y1, bounds=([0.1,0.1,peak1_disp[i]-50,0.5,0.5,-100,-1000],[1000000,1000000,peak1_disp[i]+50,30,30,100,5000]))
    I = f_psVgt(x1,popt[0],popt[1],popt[2],popt[3],popt[4],0,0)
    #cons.write('Optimized parameters:')
    #cons.write(popt)
    perr = np.sqrt(np.diag(pcov))
    #cons.write('Standard deviation:')
    #cons.write(perr)  
    p2.plot(x1, y1, pen=(255,0,0))
    p2.plot(x1, f_psVgt(x1, *popt), pen=(0,0,255))
    message =''
    for j in range(0,5):
        if abs(perr[j]/popt[j]) > 0.1:
            message = 'BAD FIT !!!'
    # fitVal.setText('A_gau:'+str(popt[0])+'+-'+str(perr[0])+'\nA_lor:'+str(popt[1])+'+-'+str(perr[1])+'\nCen:'+str(popt[2])+'+-'+str(perr[2])+'\nWidth_g:'+str(popt[3])+'+-'+str(perr[3])+'\nWidth_lor:'+str(popt[4])+'+-'+str(perr[4])+'\nBG_x:'+str(popt[5])+'+-'+str(perr[5])+'\nBG_0:'+str(popt[6])+'+-'+str(perr[6]))
    fitVal.setText(message+'\n\nA_gau: %.2f +- %.2f \nA_lor:  %.2f +- %.2f \nCen:  %.1f +- %.1f \nWidth_g:  %.3f +- %.3f \nWidth_lor:  %.3f +- %.3f \nBG_x:  %.5f +- %.5f  \nBG_0:  %.2f +- %.2f' % (popt[0],perr[0],popt[1],perr[1],popt[2],perr[2],popt[3],perr[3],popt[4],perr[4],popt[5],perr[5],popt[6],perr[6]))
    
def update_linePlot(y_index):
    global q, int_data, p2
    p2.clear()
    #p2.plot(q,int_data[y_index],pen=(255,0,0))
    #p2.plot(int_data[y_index],pen=(255,0,0))
    fit_PseudoVoigt_single(y_index)
    
def slider_update():
    global p1, p3, p4, mu, peak1_I, img1, a1, a3, a4
    p1.removeItem(a1)
    p3.removeItem(a3)
    p4.removeItem(a4)
    p5.removeItem(a51)
    p5.removeItem(a52)
    #sliderVal.setText(str(slider.value())+'/'+str(len(int_data))+' '+str(energy[slider.value()])+'eV')
    sliderVal.setText('%0.2f eV Range:[%0.2f, %0.2f]' %(energy[slider.value()], energy[0], energy[int(len(energy)/2)]))
    
    a1.setPos(0,slider.value())
    a3.setPos(energy[slider.value()],mu[slider.value()])
    if peak1_I:
        a4.setPos(energy[slider.value()],peak1_I[slider.value()])
    a51.setPos(energy[slider.value()],np.mean(I0, axis=0)[slider.value()])
    a52.setPos(energy[slider.value()],np.mean(I1, axis=0)[slider.value()])
    
    p1.addItem(a1)
    p3.addItem(a3)
    p4.addItem(a4)
    p5.addItem(a51)
    p5.addItem(a52)
    update_linePlot(slider.value())
    # win.show()

def DAFS_corrections():
    global peak1_I, mu
    p4.clear()
    peak1_I_abs = peak1_I
    if absCheckBox.checkState():
        peak1_I_abs = peak1_I*np.exp(mu)
        peak1_I_SG = savgol_filter(peak1_I_abs, savgolWin.value(), savgolPol.value())
        p4.plot(energy, peak1_I_abs, pen=(0,0,255))
    else:
        peak1_I_SG = savgol_filter(peak1_I, savgolWin.value(), savgolPol.value())
        p4.plot(energy, peak1_I, pen=(0,0,0))
        
    if savgolCheckBox.checkState():
        p4.plot(energy, peak1_I_SG, pen=(255,0,0))
    else:
        p4.plot(energy, peak1_I_abs, pen=(0,0,0))
    
    p4.show()

loadResButton.clicked.connect(load_resultData)
calcButton.clicked.connect(fit_PseudoVoigt_all)
slider.valueChanged.connect(slider_update)
logCheckBox.stateChanged.connect(log_resultData)
normCheckBox.stateChanged.connect(load_resultData)
normCheckBox2.stateChanged.connect(load_resultData)
absCheckBox.stateChanged.connect(DAFS_corrections)
savgolCheckBox.stateChanged.connect(DAFS_corrections)


#saveCurveButton.clicked.connect(update_linePlot(slider.value(),True,False))
#clrCurveButton.clicked.connect(update_linePlot(slider.value(),False,True))

win.show()

## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()