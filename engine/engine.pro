TEMPLATE = subdirs
CONFIG  += ordered
SUBDIRS += audio
SUBDIRS += src
SUBDIRS += analyzer
!android:!ios {
  SUBDIRS += test
}
