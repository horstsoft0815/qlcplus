TEMPLATE = subdirs
CONFIG  += ordered
SUBDIRS += audio
SUBDIRS += src
SUBDIRS +=
!android:!ios {
  SUBDIRS += test
}
