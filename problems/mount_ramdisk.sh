#!/bin/sh

TARGET=./ramdisk

sudo mount -t ramfs ramfs $TARGET
sudo chmod o+w $TARGET

# this may swap:
#sudo mount -t tmpfs -o size=200M none $TARGET user
