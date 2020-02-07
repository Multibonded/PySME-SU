import sys
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog
from PyQt5.QtGui import QIcon
import PyQt5.QtWidgets as qt
import matplotlib.pyplot as plt
from astropy.visualization import astropy_mpl_style


from PyQt5.QtWidgets import *
print(QStyleFactory.keys())

app = QApplication
app.setStyle('Fusion')

from astropy.io import fits
import numpy as np

# GUI initiation
class App(QWidget):
    # Sets geometry of main body
    def __init__(self):
        # need super else error
        super().__init__()
        self.title = 'Flux and wavelength plotter'
        self.left = 500
        self.top = 500
        self.width = 440
        self.height = 280
        self.layout0 = qt.QGridLayout()
        self.layout = qt.QGridLayout()
        self.group_box_settings = QGroupBox(self)
        self.group_box_settings.setTitle("Select a file to enable these many magic tricks")

        self.initUI()

    # Creates parts inside the main UI such as buttons that connect to functions when pushed
    def initUI(self):

        # Set title of window and its geometry to previously defined values
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)



        # File explorer to select a fits file to later plot. Can change to other files later if they are ever plottable
        self.filechoose = qt.QPushButton('Choose fits file to plot from')
        self.layout0.addWidget(self.filechoose, 0, 0, 1, 2)
        self.filechoose.clicked.connect(self.saveFileDialog)

        # Plot the chosen file with the file address storedfrom saveFileDialog function
        self.plotbutton = qt.QPushButton('Plot Wavelength against Flux')
        self.layout.addWidget(self.plotbutton, 1, 0, 1, 2)
        self.plotbutton.setEnabled(False)
        self.plotbutton.clicked.connect(self.PlotWF)

        # Prints all headers when pressed
        self.prnthdrsbutton = qt.QPushButton('View all headers')
        self.layout.addWidget(self.prnthdrsbutton, 2, 0, 1, 2)
        self.prnthdrsbutton.setEnabled(False)
        self.prnthdrsbutton.clicked.connect(self.prnthdrs)

        # Box to input text to search for in the headers and comments and values
        self.searchhdrsbox = qt.QLineEdit('Search for a header/comment')
        self.layout.addWidget(self.searchhdrsbox, 3, 0, 1, 1)
        self.searchhdrsbox.setEnabled(False)
        # Connects pressing enter with the search button
        self.searchhdrsbox.returnPressed.connect(self.searchhdrs)
        self.searchhdrsbox.mousePressEvent = lambda _: self.searchhdrsbox.selectAll()

        # Actually searches for the hdrs box text
        self.searchhdrsbttn = qt.QPushButton("Search")
        self.layout.addWidget(self.searchhdrsbttn, 3, 1)
        self.searchhdrsbttn.setEnabled(False)
        self.searchhdrsbttn.clicked.connect(self.searchhdrs)


        self.group_box_settings.setLayout(self.layout)
        self.layout0.addWidget(self.group_box_settings)
        self.setLayout(self.layout0)

        self.show()



        # Function for choosing file path to plot later, also enables the plot button.
    def saveFileDialog(self):
            # Parent of self, no directory or name, but restricted to fits files
            self.filePath = QFileDialog.getOpenFileName(self, '', '', '*.fits')[0]
            print(self.filePath)
            if self.filePath:
                # In case I need the name of the file itself.
                self.fileName = self.filePath.split("/")[-1]
                self.image_data = fits.open(self.filePath, ext=0)

                print(self.fileName)
                self.plotbutton.setEnabled(True)
                self.prnthdrsbutton.setEnabled(True)
                self.searchhdrsbox.setEnabled(True)
                self.searchhdrsbttn.setEnabled(True)


    # Function to plot the chosen fits file.
    def PlotWF(self):

        # the .data values of the file which is the flux in this case
        flux=(self.image_data[0].data)

        # Takes the value of CRVAL1 which is the initial angstrom wavelength
        starter =  self.image_data[0].header['CRVAL1']
        # The step in wavelength per data point taken from CDELT1
        steps = self.image_data[0].header['CDELT1']
        # Making a wavelength list with the interval of steps
        wlist = starter + (steps*np.arange(self.image_data[0].header['NAXIS1']))


        # Plot the figure with thinner lines and show it
        plt.figure()
        plt.plot(wlist,flux,linewidth=0.1)
        plt.xlabel( 'Wavelength (A)')
        plt.title(str(self.fileName))
        plt.ylabel( "Flux (Relative)")
        plt.show()

    # Prints all headers found in the file
    def prnthdrs(self):
        print (self.image_data.info())
        print(repr(self.image_data[0].header))

    # Search the headers for the text in the search box and prints entire line where its found
    def searchhdrs(self):
        print("\nSearching for '" + self.searchhdrsbox.text() + "'\n")
        for line in repr(self.image_data[0].header).split("\n"):
            if (self.searchhdrsbox.text().lower()) in (line.lower()):
                print(line)

        print("\nSearch finished\n")

class ButtonGroupBox(QWidget):

    def __init__(self, parent=None):
        super(ButtonGroupBox, self).__init__(parent=parent)

        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0,24,0,0)
        self.groupBox = QGroupBox(self)
        self.button = QPushButton("FOO", parent=self)
        self.layout.addWidget(self.groupBox)

        self.button.move(0, -4)



def main():
    app = QApplication(sys.argv)
    main = App()
    main.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()