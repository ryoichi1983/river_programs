#!/bin/sh

# ｺﾝﾊﾟｲﾙ等を行い,ﾌﾟﾛｸﾞﾗﾑが動くようにする.

make clean
chmod +rw *
dos2unix *
make
chmod +x 001FitTank.sh

