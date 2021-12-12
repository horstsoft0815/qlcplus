TEMPLATE = subdirs

!android:!ios {

 SUBDIRS              += dmxusb
 SUBDIRS              += peperoni
 SUBDIRS              += udmx
 SUBDIRS              += midi
 unix {
   system(pkg-config --exists libola) {
     system(pkg-config --exists libolaserver) {
       SUBDIRS        += ola
     }
   }
 }
 !macx:!win32:SUBDIRS += dmx4linux
 SUBDIRS              += velleman
 SUBDIRS              += enttecwing
 SUBDIRS              += hid
 !macx:!win32:SUBDIRS += spi

 greaterThan(QT_MAJOR_VERSION, 4) {
    SUBDIRS              += os2l
    #!macx:!win32:SUBDIRS += uart
    #!macx:!win32:SUBDIRS += gpio
 }
}

SUBDIRS              += artnet
SUBDIRS              += E1.31
SUBDIRS              += loopback
SUBDIRS              += osc
SUBDIRS              += qm-dsp
SUBDIRS              += samplerate
SUBDIRS              += btrack
#SUBDIRS              += btrack/modules-and-plug-ins/vamp-plugin/btrack_vamp.pro
