#!/usr/bin/vim -s
:%s/updu(0,0,b)/updu[0][0][b]/g
:%s/updu(1,0,b)/updu[1][0][b]/g
:%s/updu(2,0,b)/updu[2][0][b]/g

:%s/updu(0,1,b)/updu[0][1][b]/g
:%s/updu(1,1,b)/updu[1][1][b]/g
:%s/updu(2,1,b)/updu[2][1][b]/g

:%s/updu(0,2,b)/updu[0][2][b]/g
:%s/updu(1,2,b)/updu[1][2][b]/g
:%s/updu(2,2,b)/updu[2][2][b]/g
