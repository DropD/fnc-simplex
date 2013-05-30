#!/bin/sh

TARGET=./ramdisk

sudo mount -t ramfs ramfs $TARGET

# this may swap:
#sudo mount -t tmpfs -o size=200M none $TARGET