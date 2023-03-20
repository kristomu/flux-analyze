#!/bin/sh
# Use find to handle extremely large numbers of files.
find check_sectors/ -iname "t*" -print0 |xargs -0 rm
