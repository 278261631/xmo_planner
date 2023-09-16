(function(){"use strict";try{if(typeof document<"u"){var a=document.createElement("style");a.appendChild(document.createTextNode('.aladin-container{position:relative;border:1px solid #ddd;height:100%;overscroll-behavior-x:none}.aladin-imageCanvas,.aladin-gridCanvas{position:absolute;z-index:1;left:0;top:0}.aladin-catalogCanvas{position:absolute;z-index:2;left:0;top:0}.aladin-logo-container{position:absolute;bottom:2px;right:5px;z-index:20;min-width:32px;max-width:90px}.aladin-logo-small{padding:100% 0 0;background:url(data:image/gif;base64,R0lGODlhIAAgAJEAAJIsLdEwJAdMmP///yH5BAkAAAMALAAAAAAgACAAAAjMAAcIHEiwoMGDCBMqXMiwocOHECMaFCCxYkKKAAoK2MiRo0UBAEKKFOkxYUaCIEMSHBlyo0OQCke6HHDyJEWBKgcG2MlzoEyFMAXyHNqTZsubNFGeHLDT4FCcLREGZUqwaFGRUk82FfqUaQCoSH0OCLqVqlCuX42u9Kl1a1qzXnGGVaozLdG6cpMWxOrVblm4AOYOTNn2L1efYZdu5Eu0cV6cE0fW7QqV4WK+CAMLPnhZMtvAEDmy/CkWMtCOHVFaXC2VtevXsGPLZhgQADs=);background-size:100%;background-repeat:no-repeat}.aladin-logo-large{padding:58.45% 0 0;background:url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAI4AAABTCAMAAAB+g8/LAAACx1BMVEVMaXGBYYSWk5i7ur1fGUW0Fzbi4OP////Qz9K2s7f////qyseffX7TxczMytBXU1ndrahOWXi0o7RaH0v///+1GjfYkY29srb///+1GTe0Fzajn6RgFkFdHkni3+GLV3PU0dXMubr6+vpmIktUJVKiGDqGcX7p5ujLwMJgFkFgFkFNOWnp1tZaHUi0FzaEZohkX2VVKVXUwcvy8vI4U4tQMWBXIk+NGT9ZIEx+Wn5vF0EUYqF3c3lgFkL5+PkUYqH///////9lFkG0FzYUYqFeNF/BwMs2WpP6+vrBv8JSJ1TNy85TJlO0FzaJhYsUYqF5GEEUYqF2Zo60FzazFza0FzYUYqGWdIsrWpWTGj6jGDp3Kk58Y4S0FzZgFkFXIU2OiY+vmqVhGENlGEJqQ2z///9SKFJTJlP///9pF0GOjpd0Ol6rFzi9sbm0Fza0FzYUYqGXmLp3TXJmHkhLSXy/jJBVK1ivrLDu7e7l5OYLCw6AYYRpFkGCIUYVYqGAZoqJfofez9hZPGtcW4phFkIUYqGVbG1BToTFw8ZqZGr4+PmIGkAWYqD6+vpaHUoUYqGEZoh5ZH2ceYAbGyCmFzmgGjsUYqGAYIOuiJJ3SW1PZJlNM0OliJ+MQF5uF0Gcmp8kXZpSKFWEZojDwcXq1tQzVY9pN2CyFzbZlZFHbKOZgpWjnaRlMlsUYqGHGD9FRElaHUiZfpfW1dddW2HMtsJ3k8NTJlPDT1WlMElcGkY6UYjMa2tDSH3IpKOEZoiFTWqni54DAwQsLDGsqa3Pu8cUFBnEtr8gHyU4Nz3cwsMKDA/GV1tGRUtCKjDczM7NfXzMvcza1Nv///+9PUmhfZRxY2y2KT/15eLo4ud5fKXCXmTnu7ekZ3pgFkFTJlOEZoiUGT5aHkp8GEBzF0G0FzadGDtKQnNeJ1JqJk5fGEReGkaDGT8UYqGlSw8iAAAAwXRSTlMA87vu8R/SwN6iQP7+/vf9/J75s4DT/v0gokr33vzj++7+9/Hz8/3u1tFw9P4f5nP9cvl0/vb+/vL79HH9++WPMFA7s1r++vRhscXEiWT9PvLQ+Ffzih/9/vb+9z3Enn7N/cWI/RDWPND+9/38gTx6uPj5/fn+/efauu7k8fnl0+ro/f33wvj7meDU2PeaZquWH9jJ1O0QrPfC0vXo+uHj+J7ETZvkpfzI+6e44qCorUr22cpX3xDd9VdUvtb6V9z+sGF5dwAACP1JREFUeF7s011r01AcBvATON8gFCgkV+2AFrKSm5MGCEKlDIqCgEgpXYUaOkanQLrtpupgCxTY8A3EDWToYBNlgFeiIOIX+f/z0pe96IcwSZtRxY0ByXaT3204nIfnPCHXLJFIJBKJgoe8LLyp/+fbPXJ16mvW3k7XsjiOs3xGd+1FoVAn12Hh1g7HqcYqMsdxGAZ0K8B15avOUkGPQymFvm0Plb6InrKOuqEbqoHVd1vPSfxk+fvT/VZRpBQ0aoLPtRW7VptRKD0VGTKcmNva/0biJPmVjDZUtXN8egKBXIM3IeC64NEohHlGvV6WxOcTj4hHhmq015dHyASh0ciXSKjUhAka5in21AMSi0ev3v7UEfEEjM5Rtbd+mPssSeQfz8JEIgZoR7VIHB6ubFvj4WqQ4xvnTqIkgE+j6KPQiSHOe54vlx0Krj38BYJ08bp27UUAcZyHQibiOJIsV9DXV4a1mrKYk8jFSndn+qCJwXuJZmYt2mKy6HvyemlJ8Zd7iSO3Bx8ANKCITDONQpTVtNCzam2vfHVBOK+OvLek/FRpmy4ABWBIob0X5TsF1Th6FY/NHC9NN5BOzadvzg5m06ldmGiSiQYAOCYwBpmNHyQaX+QW+ljbPDjkH5CJheCnnx+MDZU7j+FMcyqOSDU0Ye5jNL1UshhwaNvwo4SK4mYqNQjZGvzl/lkck1GKsPz7xiUu+0Nq2b+2VYVx/NDZJTYmnV2TpuvMsiJNhbSUZmMwSpssENJl7XSmrrDNpkpn3dqO4eraoqXFMmddBWcVncImDpgOMKiiImJu3t+Wl9a54UiccOxA8keY+5xzc25ugiTx+9s5fHL55D7nPM9dk5FY6NpO1wVgJ8g0pVIpv793mWLP31JEeiMKiCa5yeu8CRIeP8STySzLIMv5VSrl+e1YLne0Ap3BMMcnNE/XdV5Ybyer+lcOZyGeIsyKn+AxSDR8qcVwq9X6Lj+sDuwlm8FMJsiJ4o2fSX9fyeeXuY2D6MrpvDz1KEtylmIG/uh2Y6ZDlOomGxBaxx86CzovybniRG12VEEMUaCXLGV03svSPPaMXsBG8jKCDssHc3aE1BgLOj9OCzoshoYKdExxYL3zpTpuODZbo6+f7hKw0A5e5sBDqQ63MGcfwkxnHZXqeL+pQEd7kbpLdY5kwebt0f1HeGwbwYy8zsGMC7Ain9UfmE5va32pDqfXVuCjCwB73Vys0wUy+0f3fV6EeWLqkRn0U13QR9MTEOql4HXI5nZE304Ilo2E6KmkWnYCh9eKdMhI2LpxwU2xaYp10lZsdWKsbj138klVD/X55Q+Mnc/mOyC0bKLjvf3c4sBJB7mX8ekKdCb0rFpMh7ThrcPCNJhRK9kVrG/txkKGkMvHQe48wOpdu1dop6Q6j6N8Glxs8R9pgNAyXDSLdIJZyE4B+zkWS4QE7Fw33oyRYKxGyEWLYVTXmz/5jn+kGY0FRQYT8kp0tJPNfDb6AI6bpDrURtt/U6PRzArYTX5IaXZo+NzDGI+g99NE5/ivu5ebIbKxv1rEBhXpmL6F0yYn1YrqpDpjFHsHsCaKJUR9JwI66Dp5cY2fHaL3SZ75p3qd1QV4yLSDlkEr0mE2XcYQYF9RbHyzSMeaR66SpnS6GcmFrvzIVq2OthMgn9YyTP6cSawj2LhPJGCnrYAlxTrOeoROXSKH52umc2FfVTqsCFE9QgagAw6RztNuavNG8i7s5DE9wSIiHesuNNONP/ZKdFS5RXm1Oqtwo8KDhbGun0DIRXUKNlNGKab8HXRo8x5xYkyP8m1LQWcAVauj1QEz/AVC5jOkDHbk7mAzi9hsklr1ibAk04GBOksb4by2y8bRn1elw2rFqWACwLwOda6/WqTjXpnCyR6GGQAL7FWfuspuFk7aomRK9L+40lKzzhwUIQBNfzAOvOpgRqxzaOVvjCMi7HJc6N91gs7DE+M+OrWW9mSequ3tsFo19svymWwjFdlT0OF3dRGFIpkog1kEnZag0hfmSO4YX9u6UrOOqYcrSWic6LB4H5TDHENwdooSMB6/AfepNh2olTTpEh1jOUyJS3QCCU/uygCqUQfmeGmGz0p0wvfLYjGpTih9/ti1F1CtOvCVU5qwR/KZd7etLDbbIcHaz+euIVS7jiPAlYsKziiLr688tsSwhU877tu+XDyK/ofOxIZMHH3KD4m0D6q2QVpINu4p8lHyiQCRUCh6lYb2tUkZRJdI+5v+fCs38BGCyGgQaofHqC7DtrD4tx07aGkbDAM4/hTmB5gFhqAILAFs0SHYpqaMwkwRhtBWtmp0FobFURqw1uJlaQdO6SVMB0zZmNCeelLmbd1p32CXIjj2BNNkZUnyIZa0tKlujAFtveR3ed/b++fhvbwv/JcvDVFDmaSQg7YzSrkhile6MjW3OQQt4Ekkxp/PhsPJmRgDvZQp3mdlXVE4Bdo8tP36pqI0z/MP8d1T6FIdVWeXxEDW9TICPRUXfFwFzRzliZ0T/UnV63XqyhqL5Y77EXR58D5dW/KryUXXIfTY6TzBss2cNTsHdVlOIVIcRSPi3vq1lmNXdrx2guF548NbgJ4PR02lsG7mjEDHKCJP0/wen5hITEK3Y5crvY1oxRRC0HMHMyparudA1T0x0SmxTbqzaTTtzhvCaRx6blLwYTtnCv5paHPkbNSKGcuVDCF4BH1QXg50cuzx/GlzZO3iG5nO1jBcNIxCEPpjoyFhE0WSCgd/88IzZ/26kT++tq6MEItAv2yI2u4YoqZpiKR+8x+9ulB+TIiSTHKsjL+aVybGHEH/lEXMhRElUULUFZ1f94DlzfT0gntjJ5kVTX5JRZ0lKyclI8NAX00TGiKqhN9cUmSF06Mpmq7L2wHRxq5UFOXzyetMKA79RgQQ0TycCEgqpnRdJ/NsXkaU8kvnH4fvnSe9Oe9qfnXZ2I/DAHwq5cY0QrT4Ec0d4feLor5y8X14a+vycnExFotlQgwMSkQo+cRWD2EuLTve3LIh7L86fAaDFr/rbRgzXsuOz+fzFnNFo3AQZODWMJmCYdsPReDWMXEm2NTd4nA4HA6H4zc5mbo+QO8AVQAAAABJRU5ErkJggg==);background-size:100%;background-repeat:no-repeat}.aladin-col{float:left;width:45%;margin-right:5%}.aladin-row:after{content:"";display:table;clear:both}.aladin-location{z-index:20;position:absolute;top:0px;padding:2px 4px;background-color:#ffffff80;font-size:11px}.aladin-projection-select{z-index:20;position:absolute;top:0px;right:48px;padding:2px 4px;background-color:#ffffff80;font-size:11px}.aladin-location-text{font-family:monospace;font-size:13px;margin-left:4px}.aladin-clipboard:before{content:" 📋";cursor:pointer}.aladin-measurement-div{z-index:77;width:100%;position:absolute;bottom:20px;font-family:monospace;font-size:12px;max-height:20%;display:flex;flex-direction:column;flex-wrap:nowrap}.aladin-measurement-div table{-ms-overflow-style:none;scrollbar-width:none;overflow-y:scroll;overscroll-behavior-x:none;background-color:#ffffffe6;border-collapse:collapse;table-layout:fixed;white-space:nowrap;height:100%;display:block;border:2px solid black}.aladin-measurement-div .tabs button{background-color:#ffffffe6;border-radius:4px 4px 0 0;float:left;border:none;outline:none;cursor:pointer;font-family:monospace;font-size:13px;text-overflow:ellipsis;text-align:left}.aladin-measurement-div table::-webkit-scrollbar{display:none}.aladin-measurement-div table thead{position:sticky;top:0;background-color:#fff;color:#000}.aladin-measurement-div table td,.aladin-measurement-div table th{padding:.3em;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;text-align:left;color:#000}.aladin-marker-measurement{max-height:130px;overflow-y:auto;overflow-x:hidden;font-family:monospace;font-size:11px;color:#000}.aladin-marker-measurement table{table-layout:fixed;width:100%}.aladin-marker-measurement table tr td{word-wrap:break-word}.aladin-marker-measurement tr:nth-child(even){background-color:#ddd}.aladin-marker-measurement td:first-child{font-weight:700}.aladin-fov{z-index:20;position:absolute;font-size:12px;font-weight:700;font-family:monospace;bottom:0px;padding:2px;color:#321bdf;background-color:#ffffff80}.aladin-maximize{position:absolute;top:6px;right:3px;z-index:20;width:30px;height:30px;background-image:url(data:image/gif;base64,R0lGODlhFgAWAOMJAAAAAAEBAQICAgMDAwUFBAUFBQcHBwgICAsLCv///////////////////////////yH5BAEKAA8ALAAAAAAWABYAAARm8MlJabr41i0T+CCQcJz3CYMwklsSEAUhspVnHEYw0x1w558VzYQTBXm24sgjJCUQykmT9dzxWlNq9vrwILYtqe/wRc6SBqszOE6DLZ/AT00FyKNcD4wQeLdQAiB+cCFHVxkZXA8RADs=);background-repeat:no-repeat;background-position:center center}.aladin-layersControl-container{position:absolute;top:30px;left:4px;cursor:pointer;z-index:20;background:rgba(250,250,250,.8);border-radius:4px}.aladin-layersControl-container:hover{background:rgba(220,220,220,.8)}.aladin-layersControl{width:32px;height:34px;background-image:url(data:image/gif;base64,R0lGODlhGQAcAMIAAAAAADQ0NKahocvFxf///wAAAAAAAAAAACH5BAEKAAcALAAAAAAZABwAAANneLoH/hCwyaJ1dDrCuydY1gBfyYUaaZqosq0r+sKxNNP1pe98Hy2OgXBILLZGxWRSBlA6iZjgczrwWa9WIEDA7Xq/R8d3PGaSz97oFs0WYN9wiDZAr9vvYcB9v2fy/3ZqgIN0cYZYCQA7);background-repeat:no-repeat;background-position:center center}.aladin-cooGridControl-container{position:absolute;left:4px;cursor:pointer;z-index:20;background:rgba(250,250,250,.8);border-radius:2px}.aladin-cooGridControl-container:hover{background:rgba(220,220,220,.8)}.aladin-cooGridControl{width:32px;height:34px;background-image:url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABwAAAAcCAQAAADYBBcfAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAAAmJLR0QA/4ePzL8AAAAJcEhZcwAAkV8AAJFfAZ/wQlIAAAAHdElNRQfmBwYSADnbJS7lAAADc0lEQVQ4y23Ve0wWZBTH8c/LJUVJ9GUvhM0SDUtxplmRSmXmrcx1cfAyzC7TOV3TLHPaZc5c5ZzOVa6y2cywWqnV/KMsm2aWIBEkCBIm5i10oESSKaS+/QEiYue/85zz3fPbOec5T0BHCxlpggw9xQs45agSX9qm7vK06Mu8Hh633BMSHRQUUCnRVhFhk3RS7Yz/tQxb1MozUU9z7ZUtZK1Ngq4zW4kthl0JRcl2yFmTwSCVckEvheaA3tY7LFugBQi0YlO94G8JVuohxQD97VTvuAq3GGGqG91vsET/mmf9JTBslRinRan0h2PGOWKrHlL0li6oVrNNtqmywmBZ8ltuHahKkcWKvaazANZYDhIsUOOocw56Vnes02izIFFizdEs20L7NDsrgl/cKkG6PDmeN9NBHwj73IP6WCYoB0aotQAsk9cqPVWR1Yrk6YtkxSZKtkKdA9LM8qOkKI8IaQS7DJQCDtlpmgbzVeOkwwapVeqcOLm+FW9UtCVS/OZr1AtrUoRcU23U1zhdRYuT4Q5jZVsmz1wxLkjlhFLbJYAZygwwwX4zREnzsh+UK3fMn94yVBTGqvCrMho9p8w4EC9PoQrz20bxan2le0aB7i52fomIBk4ZbaXPdAG3qVPjzg5zNUaJZBAjxx4/a6ZJlv4qTENnq2z3vkovuv7icCFTmV6iDPSOKtON1RCj2VAbvGGe/foYaYpSD3tKjgI7VDnprJB4j+lvqEpPyveAJioUCYr1qiq/m972wB7yrnzlypWqdsYuK9wjDjxtD2s0y21t8z9mtYZaypAgze2GmalSv7bzFJutJcsZhVJMscdCZTbKdFWH4tynRBLoapJdGmWRJF/Ep4rNwQCrVdkgV6pObeCjdkrSz3RblKtWICmAGd4Uo9J4RxBriLC7xKlxQI1TIu4x3F4hf/lKoaXe83YAQZ8Y6bgvvNK6kgISpbtZmmt0ETHYUXl22+u0FW4SVt8iZbh9PrbDVveK7bBSYnT1jbmt3mxVhrdPCNtniddVWmeMbpfB19ptPG6wyH7h9juHgGyLbPOd0TI1+EmxA+o1ici01Ed6udt586wXaQ+2CF4swYcq9JEhTTcBF9BTSMB533tJwaUmt7egHJN1UWK3E6IlSJZqjG5KrbXhYkmuBCFklAmGCIkR0axOic22d/wC/gNPgB2kPbE+ygAAACV0RVh0ZGF0ZTpjcmVhdGUAMjAyMi0wNy0wNlQxODowMTowNSswMjowMOuc89QAAAAldEVYdGRhdGU6bW9kaWZ5ADIwMjItMDctMDZUMTg6MDA6NTcrMDI6MDCqfD8bAAAAGXRFWHRTb2Z0d2FyZQB3d3cuaW5rc2NhcGUub3Jnm+48GgAAAABJRU5ErkJggg==);background-repeat:no-repeat;background-position:center center}.aladin-layerBox{top:30px;padding:4px;max-height:90%;overflow-y:auto}.aladin-overlay-list{padding-bottom:4px}.aladin-box{display:none;z-index:30;position:absolute;left:40px;background:whitesmoke;font-size:13px;font-family:Verdana,Geneva,Tahoma,sans-serif;background:#fff;line-height:1.3;color:#222;border-radius:2px;box-shadow:0 0 6px #0003;-ms-overflow-style:none;scrollbar-width:none;overflow-y:scroll}.aladin-box::-webkit-scrollbar{display:none}.aladin-dialog{display:none;z-index:30;position:absolute;top:50%;left:50%;transform:translate(-50%,-50%);background:#eee;border-radius:4px;font-family:Verdana,Geneva,Tahoma,sans-serif;line-height:1.3;color:#222;max-width:500px;padding:.8em}.aladin-layer .aladin-input,.aladin-layer .aladin-selector{width:100%}.aladin-selector,.aladin-input{padding:4px;box-sizing:border-box;margin:0}.aladin-layerSelection{max-width:150px}canvas{image-rendering:optimizeSpeed;image-rendering:-moz-crisp-edges;image-rendering:-webkit-optimize-contrast;image-rendering:-o-crisp-edges;image-rendering:pixelated;-ms-interpolation-mode:nearest-neighbor}.aladin-box-title{font-size:16px;font-family:Verdana,Geneva,Tahoma,sans-serif;line-height:2;font-weight:700;cursor:pointer}.aladin-box-content{padding:10px}.aladin-text{background-color:#bbb}.aladin-gotoControl-container{position:absolute;left:4px;cursor:pointer;z-index:20;background:rgba(250,250,250,.8);border-radius:2px}.aladin-gotoControl-container:hover{background:rgba(220,220,220,.8)}.aladin-gotoControl{width:32px;height:34px;background-image:url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAADgdz34AAAAGXRFWHRTb2Z0d2FyZQBBZG9iZSBJbWFnZVJlYWR5ccllPAAAASZJREFUeNpiYMAOpIG4GogPAPETIP4PpQ9AxaUZyARsQNwFxN+ghuLC36Dq2EgxXBSITxIwGB2fhOojyuXohl8C4hCk4JCG8i9hsYSgT9qRNPwE4hY8mtig8n+Q9LTjM1waaihMcQ6RQZqD5ig5XAobkRQeJjFRHEbS20iMoigSLYgixnHPkRSRmr6lkfQ+x6UIOfzZyMg3yPFAfx8wAfFNJL49iRYgq7+JS1EFlVJRBS5FcmjxkEak4WnE5gMGaMGFrBhUYjLjUMsMlUd21CJyyqLz0LJHAqpGAso/j6XQ+0NMHiKnNEW3JJyQJZzQ4PpJwLCf0GD5Q44l6DXac6R0jl6jhVBiCbEAlyUh9LDEm9aWPGegMkC35AkDDYAf1OUg7AoQYADEj7juEX/HNAAAAABJRU5ErkJggg==);background-repeat:no-repeat;background-position:center center}.aladin-gotoBox{top:68px;padding:4px}.aladin-target-form{padding:5px}.aladin-simbadPointerControl-container{position:absolute;left:4px;cursor:pointer;z-index:20;background:rgba(250,250,250,.8);border-radius:2px}.aladin-simbadPointerControl-container:hover{background:rgba(220,220,220,.8)}.aladin-simbadPointerControl{width:32px;height:34px;background-image:url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAYAAADgdz34AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEIDgUy0LKZuQAABc5JREFUSEuFln9sVeUZxz/P+97TbmtL0pZrO5wtvT/bkolBQ5gJrh1FO9RNTYSEYEoiywjiHBlbssXNP4Y2ziZLlBJ1GBFnWF2CPya0TgxsMytmY8Fk4/bec3sv1AXaXmYyaBvovec8++Oe2xUw25uc5LzvOe/3+fV9fkhPT1eI/y4BdGEjIs03N3kHXz3kiwix5W0/ERF187l+VZW+RzfLhU8njapec28xjiwSsHCoqqy5a7X/9FPP+R3xRNjzvGagRlW3qyrGmBeBWRuyk6lMpvDm8Ovy0sB+KyI3KhkIkDKu4lQ5jAwfK3XEE02lUmmNqm4AuoE41y4XOC4iR23Inhxz3ane3p5QsVgkELQgwKmA1y2p1bcOv+clo7E7fN9/DNga/Pe3ADAW7LOUBa4K9geMMYOZ3PhfH3jwPnv50mWpWGMjkeUGwKly+N27I14iGv2G+joI9AInRORZERmMJ1te/Oxf/24UkVRrZNkPLl2aGRXkE6AGeFBVVzc2NKT/PPpxPpmMG9/3BVDp6ekKqSoffviHUiISvcP3/f3ASmCfDdlfprPZbKAl8bbIkwBuPrencpaMxmKe5+0CdgCnjTHbMrnxU+vWfT0kIhiANXet9jviiabALSuBfdVfqP5pOpvNqmoFKKqqSVVNJqOxCMBvR34t6fFstqq66klgH3Cb7/s72+Pxpke2b/YAZP36bueDD44X422Rb6vq28AJG7LfqWjeEU/cVCwWnwPCwSPAFHDRcZzdKTdTWGTJr4AuEXnAzefeWb++2zHNNzd5HfFEOGALIjJU0TwAHwImjDW7gT8Bf7TW/gg4WywWhzriifC5uTFJj2ezIjIEoKobOpPJcLg57JuDrx7yA553U2bLaCCIQPOPnKrQ85nx8TNAQYRCejz7D6fK2Qt8VCwWB1q/1F7h/WiA0VUqlpreODjkm4BONZRp58YSt/w9MDkKhI01h1IZt7Dh3rtt3ZLaF2pqa57fcO/dNpXJFIw1vwHClZg03XLTGcp0TgC1IqISpH9CVfuAU8BhETGqmgTaKbulULek9oWlTUuvAExfmK6enZl7QkSWqupaYAxII/goDwG3A68BGSMiWmEKcE0WUg4oIkEmLvrnc+6pIFK5A6gIatx8rj+oLQDZlraWX7j53B5r7VPAlLX2FTef71/WsuzK9IXp6sJkofrLX2m+6uZz/caaV4CLNmR/5uZzT7e2tTxL2UUYY1528/n+SiWcDT7EJ85NrABIj2dzwEXP8zZ2JBLho0d+783OzD0xOzP3veHhY6XORCLse/4mYDqdzeYBzuUnVhDEErisqmL6Ht0sNmQngePAKny+BuUkchxnN7C2OF98PBGNdgY+Dyei0RXz88XHgbWO4/xwkavupFyfTlhrp7Zu2yLmwqeTJpXJFETkKICqbkxGY7GHe7doys0UHMfZBLT6nj8QBHSt7/kDQKvjOJtSbmZaREhGYzFV3QggIkdSbqZwfuK8Maoqbw6/LjZkTwIHgG7P83Z1JpP1ACOfvHsxezbfZ63dScAWG7I7smfzfWcy6WmAFe3t9Z7nfR/oAg7YkD05MPiMqKoYQF8a2G/HXHfKGDMInAZ2zF+d35OMxmKVJApikhaRVMXnFc2vXrn6c+AxysVu75jrTg0ffj8EqI1ElhsRIZGImY//cuqfjQ0NaZTbgPtV9dbG+oYvNtbXFztvTX42c2n2TkFk1aqVoxbz1cb6+od83/8xsBE4LSK73Hxu9J7edaFKXBY6mqpq3ZI6ffut97xEJHq77/s7+fyGI8H79Q1nbyY3fur+b30zNDc7t5BPCx0Nyr3YcRxGRo6V2uPxJq/k/b+WeUJEjlRa5j2960KlYmlxsopcN1VAYM0j2zd7Wx/+rnYmk+FSsdQM1Pi+vx1QY8zLwGVr7VTKzRQGBp+R4cPvh0TkhsniegE3jC3h5rD/xsEhX0S0XLdQN5/vV1XZum2LnJ84/z/Hlv8AT3rLVQaE2oEAAAAASUVORK5CYII=);background-repeat:no-repeat;background-position:center center}.aladin-shareControl-container{position:absolute;bottom:30px;left:4px;cursor:pointer;z-index:20;background:rgba(250,250,250,.8);border-radius:2px}.aladin-shareControl-container:hover{background:rgba(220,220,220,.8)}.aladin-shareControl{width:32px;height:34px;background-image:url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAAYCAMAAADXqc3KAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAA/1BMVEUAAAAzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzMzP////6PkQMAAAAU3RSTlMAAAI2l9TiyXwcBmrj/cEvePX69MD+Vzl06Px5sMgDZPsbOwUfKFGAHVjSEkBgXTcMllbG8rbhdtcm26Ko5/M8F3clPqN1JK+usRpVQdNQBGN9OkYs7f8AAAABYktHRFTkA4ilAAAACXBIWXMAARCQAAEQkAGJrNK4AAABMUlEQVQoz13SeVuCQBAGcAYV5TREEpUQr9TMK9G8Nc/usvn+36VdRFDmz9/uwzuzA8NcFrCRaIyLJ3iAaxdECUnJShKu/UZFTGkyYpoPu67cZows5hKXnif3TfJ5MHS8C9wqcDaaRRILJQ254KBcqd7X6tShYePDmaH52Gp3uq7zacTe+fpTnwQ7gwxJKJqIw4iX++ygW6Ox66oIQf+T6WyOi+WK+osAvptri91sd/tDyF/3u+2GtdZvIT+slgucz6YT3xkQqb/DeHSKd3YnZzpD1wEyA0dHtf9R9h62R5/Snavbabc+m/4mvtBuuF6vfVcrP8EeONRK1MlcNvdrBQdx1A2SwNO58sLF4o45zBp/DUX3+/SKV8iOtRSGnYGkItP+JVEI/RnAH+NSLBphr/0f6ns3OQ1Zz7kAAAAldEVYdGRhdGU6Y3JlYXRlADIwMTMtMTEtMThUMTQ6MzM6NTUrMDE6MDCINg3ZAAAAJXRFWHRkYXRlOm1vZGlmeQAyMDEzLTExLTE4VDE0OjMzOjA0KzAxOjAwF/ywtQAAABl0RVh0U29mdHdhcmUAd3d3Lmlua3NjYXBlLm9yZ5vuPBoAAAAASUVORK5CYII=);background-repeat:no-repeat;background-position:center center}.aladin-shareBox{bottom:30px;padding:4px}.aladin-shareInput{width:300px;margin-top:10px;display:block}.aladin-target-form input{width:140px;margin-left:5px}.aladin-cb-list ul{margin:0;padding-left:4px;list-style:none;display:flex;flex-direction:column;align-items:stretch}.aladin-cb-list ul li{flex:1 1 100%;text-align:left}.aladin-cb-list label{display:inline-block}.aladin-cb-list input[type=checkbox]{margin:3px}.aladin-indicatorBtn{color:#605f61!important;border:1px solid #AEAEAE!important;border-radius:3px;font-size:13px!important;background:#fff!important;margin-right:.2em!important;font-family:monospace;cursor:pointer}.aladin-indicatorBtn:hover{color:#201f21!important}.aladin-closeBtn{float:right;margin-top:0 0 2px 0;cursor:pointer;color:#605f61;border:1px solid #AEAEAE;border-radius:3px;background:#fff;font-size:15px;font-weight:700;display:inline-block;line-height:0px;padding:8px 2px}.aladin-closeBtn:hover{color:#201f21}.aladin-label{font-weight:700;padding-bottom:4px;font-family:Verdana,Geneva,Tahoma,sans-serif}.aladin-box-separator{height:0;border-top:1px solid #f2f2f2;margin:5px 0 5px -4px}.aladin-blank-separator{height:10px}.aladin-restore{position:absolute;top:6px;right:3px;z-index:20;width:30px;height:30px;background-image:url(data:image/gif;base64,R0lGODlhFgAWAOMJAAAAAAEBAQICAgMDAwUFBAUFBQcHBwgICAsLCv///////////////////////////yH5BAEKAA8ALAAAAAAWABYAAARk8MlJ60zJapsAyhsFPp1xfKHUZeVhAOPWwYD5xjIABDacXrseLtQhFAiBIVEwEOh8P9VzquRgrhgrNhudTaHdJ1NQ1SQCRgL41zIE1sSa6w0361y0eqXtW7EReCNlIgh6XYMbEQA7);background-repeat:no-repeat;background-position:center center}.aladin-fullscreen{position:fixed!important;z-index:9999;top:0;left:0;height:100%!important;width:100%!important;border:0!important;max-width:none!important;background:#fff}.aladin-zoomControl{z-index:20;position:absolute;top:50%;height:48px;right:8px;padding:0;margin:-24px 0 0;font-weight:700;font-size:18px;font-family:Verdana,Geneva,Tahoma,sans-serif}.aladin-zoomControl a{width:20px;height:20px;line-height:18px;display:block;background-color:#fafafacc;margin:1px;text-align:center;border-radius:4px;border:1px solid #aaa;text-decoration:none;color:#222}.aladin-zoomControl a:hover{background-color:#dcdcdccc}.aladin-cmSelection{width:60px;margin-right:10px}.aladin-btn{display:inline-block;padding:6px 8px;margin-bottom:0;font-size:12px;font-weight:400;text-align:center;white-space:nowrap;vertical-align:middle;cursor:pointer;border:1px solid transparent;border-radius:3px;color:#fff;background-color:#428bca;border-color:#357ebd;margin-right:4px}.aladin-layer{padding:0;margin-bottom:4px;background-color:#f2f2f2;border-radius:2px}.aladin-cancelBtn{background-color:#ca4242;border-color:#bd3935}.aladin-box-content>div{margin:10px 0 0}.aladin-btn-small{display:inline-block;border-radius:3px;margin-bottom:0;padding:0;vertical-align:middle;cursor:pointer;border:1px solid transparent;color:#fff;font-size:14px;background-color:#428bca;border-color:#357ebd;margin-left:2px;min-width:1.5em;min-height:1.5em}.aladin-button:hover{color:#fff;background-color:#3276b1;border-color:#285e8e}.aladin-unknownObject{border:3px solid red}.aladin-popup-container{z-index:25;position:absolute;width:200px;display:none;line-height:1.3}.aladin-popup{font-family:Verdana,Geneva,Tahoma,sans-serif;font-size:13px;background:white;border:1px solid #bbb;border-radius:2px;padding:4px;top:80px;left:110px}.aladin-popup-arrow{display:block;border-color:#fff transparent transparent;border-style:solid;border-width:12px;width:0px;height:0px;margin-top:-1px;margin-left:auto;margin-right:auto}.aladin-popupTitle{font-weight:700}.aladin-options{margin-top:4px}.aladin-options .row{display:flex;align-items:center;margin-top:5px}.aladin-options .row label{text-align:center}.aladin-options .row .col-label{flex:.3}.aladin-options .row .col-input{flex:.7}.aladin-options .row .col-input select{width:100%}.aladin-layer-label{padding:0 4px;color:#ddd;border-radius:8px;cursor:pointer;user-select:none;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;vertical-align:text-top;max-width:150px}.aladin-sp-title a{text-decoration:none;color:#317d8d}.aladin-sp-content{font-size:12px}.aladin-sp-content a{text-decoration:none;color:#478ade;font-size:11px}.aladin-cuts{width:8em}.aladin-stack-icon{width:16px;height:16px;background-repeat:no-repeat;background-position:center center;display:inline-block;vertical-align:text-top}.aladin-chevron{display:inline-block;width:16px;height:16px;cursor:pointer;vertical-align:middle;background:url(data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBlbmNvZGluZz0iVVRGLTgiIHN0YW5kYWxvbmU9Im5vIj8+Cgo8c3ZnIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyIgdmlld0JveD0iMCAwIDQ0IDQ0IiBzdHJva2U9IiM0NDQiIGJhY2tncm91bmQtY29sb3I9ImJsdWUiPgoKICAgIDxyZWN0IHg9IjIiIHk9IjIiIGhlaWdodD0iNDAiIHdpZHRoPSI0MCIgZmlsbD0idHJhbnNwYXJlbnQiIHN0cm9rZS13aWR0aD0iMyIgc3Ryb2tlLWxpbmVjYXA9InJvdW5kIiBzdHJva2UtbGluZWpvaW49InJvdW5kIiAvPgogICAgPHBhdGggZD0iTTEwIDEzIEwyMiAyMSBMMzQgMTMiIGZpbGw9InRyYW5zcGFyZW50IiBzdHJva2Utd2lkdGg9IjYiIHN0cm9rZS1saW5lY2FwPSJyb3VuZCIgc3Ryb2tlLWxpbmVqb2luPSJyb3VuZCIvPiAKICAgIDxwYXRoIGQ9Ik05IDI1IEwyMiAzMyBMMzUgMjUiIGZpbGw9InRyYW5zcGFyZW50IiBzdHJva2Utd2lkdGg9IjYiIHN0cm9rZS1saW5lY2FwPSJyb3VuZCIgc3Ryb2tlLWxpbmVqb2luPSJyb3VuZCIvPiAKCjwvc3ZnPgo=) no-repeat}.aladin-chevron-down{transform:rotate(0)}.aladin-chevron-left{transform:rotate(90deg)}.aladin-chevron-up{transform:rotate(180deg)}.aladin-chevron-right{transform:rotate(270deg)}.aladin-layer .right-triangle:before{content:"▶"}.aladin-layer .down-triangle:before{content:"▼"}.aladin-context-sub-menu,.aladin-context-menu{position:fixed;background:#fff;color:#000;z-index:9999999;width:150px;margin:0;padding:2px 0;border-radius:2px;box-shadow:0 0 6px #0003;font-size:12px;font-family:Verdana,Geneva,Tahoma,sans-serif}.aladin-context-menu .aladin-context-menu-item{height:22px;display:flex;align-items:center;padding:4px 8px;cursor:pointer;position:relative;border-bottom:1px solid #f2f2f2}.aladin-context-menu-item:hover{background:rgb(212,231,250)}.aladin-context-menu .aladin-context-menu-item span{display:block;white-space:nowrap;overflow:hidden;text-overflow:ellipsis}.aladin-context-menu-item:last-of-type{border-bottom:none}.aladin-context-menu .context-menu-item:hover{background:#f2f2f2}.aladin-context-menu .aladin-context-sub-menu{position:absolute;top:0;left:100%;display:none;width:150px}.aladin-context-menu .aladin-context-menu-item:hover>.aladin-context-sub-menu{display:block}.aladin-context-menu.left .aladin-context-sub-menu{left:0;transform:translate(-100%)}.aladin-context-menu.top .aladin-context-sub-menu{top:100%;transform:translateY(-100%)}.aladin-context-menu.left.top .aladin-context-sub-menu{transform:translate(-100%,-100%)}.aladin-sp-cursor{cursor:auto;cursor:url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QAAAAAAAD5Q7t/AAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4gEIEwUMBS20MQAAApRJREFUWEfdV7GO2zAMfXY7x1luzOFwDeB/KHCAB87eumvqmI/JWhTwF3TIrMFAv8JbgWyHG3IpOjbpcKJC02KipECHclIU+j2KfKRlAABRUyFYak3UgKiZiXUJAIt1d1ysu2PYn53DsPCucpbkRE3JAdxKTtQgWo6zJAfGGbAwLuEV7OB9/yrXIbqZ9/0+rMuhdb+RYfWmm+fged8fihzyoXXxlNcYB2KRxwwA6SzI9Crgd973B/5/u3IF8FYSw98M5D2QR86kRE3pfX+AMqKm8iEQoqYaWrfj/4bW7QgokMhCeSs5KeVrDM4K29C6Y6oE0YGEUlndocVGyg/rGZDXBQpr0kmjLiAlOOvk3vf7RJZM5Q+te2W/7coVEq+wnC1yCZYyS3B8MMYF3vBKK9Jz5PWmm8sab1eu4N9D63aaHEAsH88SzmQpnSHsHDmnWvlXIoiJ4JTviY8FRUJ8ZAiOEiJbLpcjARM1M0twpmjZOeGQrXYm1/7yMLz/eH8fnyus6XXOtmLgWBrKtdMguNKYPPystIausn9egsVDfI674KxS6003F/sX357AdOBwcADwWD9EjBKJubxYd0cml4RhP5IPrduFk04GjiZXMyZiRHKdAU3OWeB2lLNAk+tpl9IHY0zSlBo4OgsaTJpFzs+F0pxKYEXK73Nd8+3KFfWmG71qGfgSeQqPN0fvdlZrSvmUoXaJJ7EiocArZaS6tgzA+5RQe3gumUl5gPLXzy+6pEQNzDuhvFKx1ZuuEuTxliwmY/KuoHxGhzGv5QDw9PTx7senz88aLHUpJeOu8OHb17u+//4CTMm97/eFVQK5vqR8y6xOyipByjn3+0C3moUH4O++D1UX3PR9eJWzJCf6Tz5O/wDUeIfTLPlbywAAAABJRU5ErkJggg==) 15 15,auto}.autocomplete{background:white;z-index:20000;font:14px/22px -apple-system,BlinkMacSystemFont,Segoe UI,Roboto,Helvetica Neue,Arial,sans-serif;overflow:auto;box-sizing:border-box;border:1px solid rgba(50,50,50,.6)}.autocomplete *{font:inherit}.autocomplete>div{padding:0 4px}.autocomplete .group{background:#eee}.autocomplete>div:hover:not(.group),.autocomplete>div.selected{background:#81ca91;cursor:pointer}')),document.head.appendChild(a)}}catch(o){console.error("vite-plugin-css-injected-by-js",o)}})();
var wQ = typeof globalThis < "u" ? globalThis : typeof window < "u" ? window : typeof global < "u" ? global : typeof self < "u" ? self : {};
function tQ(B) {
  return B && B.__esModule && Object.prototype.hasOwnProperty.call(B, "default") ? B.default : B;
}
var sQ = { exports: {} };
/*!
 * jQuery JavaScript Library v3.6.1
 * https://jquery.com/
 *
 * Includes Sizzle.js
 * https://sizzlejs.com/
 *
 * Copyright OpenJS Foundation and other contributors
 * Released under the MIT license
 * https://jquery.org/license
 *
 * Date: 2022-08-26T17:52Z
 */
(function(B) {
  (function(A, g) {
    B.exports = A.document ? g(A, !0) : function(C) {
      if (!C.document)
        throw new Error("jQuery requires a window with a document");
      return g(C);
    };
  })(typeof window < "u" ? window : wQ, function(A, g) {
    var C = [], I = Object.getPrototypeOf, E = C.slice, o = C.flat ? function(Q) {
      return C.flat.call(Q);
    } : function(Q) {
      return C.concat.apply([], Q);
    }, t = C.push, s = C.indexOf, M = {}, N = M.toString, k = M.hasOwnProperty, U = k.toString, Y = U.call(Object), K = {}, l = function(i) {
      return typeof i == "function" && typeof i.nodeType != "number" && typeof i.item != "function";
    }, d = function(i) {
      return i != null && i === i.window;
    }, q = A.document, W = {
      type: !0,
      src: !0,
      nonce: !0,
      noModule: !0
    };
    function T(Q, i, D) {
      D = D || q;
      var w, e, c = D.createElement("script");
      if (c.text = Q, i)
        for (w in W)
          e = i[w] || i.getAttribute && i.getAttribute(w), e && c.setAttribute(w, e);
      D.head.appendChild(c).parentNode.removeChild(c);
    }
    function v(Q) {
      return Q == null ? Q + "" : typeof Q == "object" || typeof Q == "function" ? M[N.call(Q)] || "object" : typeof Q;
    }
    var _ = "3.6.1", a = function(Q, i) {
      return new a.fn.init(Q, i);
    };
    a.fn = a.prototype = {
      // The current version of jQuery being used
      jquery: _,
      constructor: a,
      // The default length of a jQuery object is 0
      length: 0,
      toArray: function() {
        return E.call(this);
      },
      // Get the Nth element in the matched element set OR
      // Get the whole matched element set as a clean array
      get: function(Q) {
        return Q == null ? E.call(this) : Q < 0 ? this[Q + this.length] : this[Q];
      },
      // Take an array of elements and push it onto the stack
      // (returning the new matched element set)
      pushStack: function(Q) {
        var i = a.merge(this.constructor(), Q);
        return i.prevObject = this, i;
      },
      // Execute a callback for every element in the matched set.
      each: function(Q) {
        return a.each(this, Q);
      },
      map: function(Q) {
        return this.pushStack(a.map(this, function(i, D) {
          return Q.call(i, D, i);
        }));
      },
      slice: function() {
        return this.pushStack(E.apply(this, arguments));
      },
      first: function() {
        return this.eq(0);
      },
      last: function() {
        return this.eq(-1);
      },
      even: function() {
        return this.pushStack(a.grep(this, function(Q, i) {
          return (i + 1) % 2;
        }));
      },
      odd: function() {
        return this.pushStack(a.grep(this, function(Q, i) {
          return i % 2;
        }));
      },
      eq: function(Q) {
        var i = this.length, D = +Q + (Q < 0 ? i : 0);
        return this.pushStack(D >= 0 && D < i ? [this[D]] : []);
      },
      end: function() {
        return this.prevObject || this.constructor();
      },
      // For internal use only.
      // Behaves like an Array's method, not like a jQuery method.
      push: t,
      sort: C.sort,
      splice: C.splice
    }, a.extend = a.fn.extend = function() {
      var Q, i, D, w, e, c, G = arguments[0] || {}, R = 1, n = arguments.length, L = !1;
      for (typeof G == "boolean" && (L = G, G = arguments[R] || {}, R++), typeof G != "object" && !l(G) && (G = {}), R === n && (G = this, R--); R < n; R++)
        if ((Q = arguments[R]) != null)
          for (i in Q)
            w = Q[i], !(i === "__proto__" || G === w) && (L && w && (a.isPlainObject(w) || (e = Array.isArray(w))) ? (D = G[i], e && !Array.isArray(D) ? c = [] : !e && !a.isPlainObject(D) ? c = {} : c = D, e = !1, G[i] = a.extend(L, c, w)) : w !== void 0 && (G[i] = w));
      return G;
    }, a.extend({
      // Unique for each copy of jQuery on the page
      expando: "jQuery" + (_ + Math.random()).replace(/\D/g, ""),
      // Assume jQuery is ready without the ready module
      isReady: !0,
      error: function(Q) {
        throw new Error(Q);
      },
      noop: function() {
      },
      isPlainObject: function(Q) {
        var i, D;
        return !Q || N.call(Q) !== "[object Object]" ? !1 : (i = I(Q), i ? (D = k.call(i, "constructor") && i.constructor, typeof D == "function" && U.call(D) === Y) : !0);
      },
      isEmptyObject: function(Q) {
        var i;
        for (i in Q)
          return !1;
        return !0;
      },
      // Evaluates a script in a provided context; falls back to the global one
      // if not specified.
      globalEval: function(Q, i, D) {
        T(Q, { nonce: i && i.nonce }, D);
      },
      each: function(Q, i) {
        var D, w = 0;
        if (iA(Q))
          for (D = Q.length; w < D && i.call(Q[w], w, Q[w]) !== !1; w++)
            ;
        else
          for (w in Q)
            if (i.call(Q[w], w, Q[w]) === !1)
              break;
        return Q;
      },
      // results is for internal usage only
      makeArray: function(Q, i) {
        var D = i || [];
        return Q != null && (iA(Object(Q)) ? a.merge(
          D,
          typeof Q == "string" ? [Q] : Q
        ) : t.call(D, Q)), D;
      },
      inArray: function(Q, i, D) {
        return i == null ? -1 : s.call(i, Q, D);
      },
      // Support: Android <=4.0 only, PhantomJS 1 only
      // push.apply(_, arraylike) throws on ancient WebKit
      merge: function(Q, i) {
        for (var D = +i.length, w = 0, e = Q.length; w < D; w++)
          Q[e++] = i[w];
        return Q.length = e, Q;
      },
      grep: function(Q, i, D) {
        for (var w, e = [], c = 0, G = Q.length, R = !D; c < G; c++)
          w = !i(Q[c], c), w !== R && e.push(Q[c]);
        return e;
      },
      // arg is for internal usage only
      map: function(Q, i, D) {
        var w, e, c = 0, G = [];
        if (iA(Q))
          for (w = Q.length; c < w; c++)
            e = i(Q[c], c, D), e != null && G.push(e);
        else
          for (c in Q)
            e = i(Q[c], c, D), e != null && G.push(e);
        return o(G);
      },
      // A global GUID counter for objects
      guid: 1,
      // jQuery.support is not used in Core but other projects attach their
      // properties to it so it needs to exist.
      support: K
    }), typeof Symbol == "function" && (a.fn[Symbol.iterator] = C[Symbol.iterator]), a.each(
      "Boolean Number String Function Array Date RegExp Object Error Symbol".split(" "),
      function(Q, i) {
        M["[object " + i + "]"] = i.toLowerCase();
      }
    );
    function iA(Q) {
      var i = !!Q && "length" in Q && Q.length, D = v(Q);
      return l(Q) || d(Q) ? !1 : D === "array" || i === 0 || typeof i == "number" && i > 0 && i - 1 in Q;
    }
    var IA = (
      /*!
       * Sizzle CSS Selector Engine v2.3.6
       * https://sizzlejs.com/
       *
       * Copyright JS Foundation and other contributors
       * Released under the MIT license
       * https://js.foundation/
       *
       * Date: 2021-02-16
       */
      function(Q) {
        var i, D, w, e, c, G, R, n, L, p, V, H, u, BA, cA, QA, TA, bA, wg, UA = "sizzle" + 1 * /* @__PURE__ */ new Date(), eA = Q.document, Eg = 0, YA = 0, fA = dI(), aI = dI(), SI = dI(), tg = dI(), fg = function(y, F) {
          return y === F && (V = !0), 0;
        }, ug = {}.hasOwnProperty, ig = [], Ug = ig.pop, hg = ig.push, Sg = ig.push, XB = ig.slice, xg = function(y, F) {
          for (var J = 0, m = y.length; J < m; J++)
            if (y[J] === F)
              return J;
          return -1;
        }, QB = "checked|selected|async|autofocus|autoplay|controls|defer|disabled|hidden|ismap|loop|multiple|open|readonly|required|scoped", lA = "[\\x20\\t\\r\\n\\f]", mg = "(?:\\\\[\\da-fA-F]{1,6}" + lA + "?|\\\\[^\\r\\n\\f]|[\\w-]|[^\0-\\x7f])+", zB = "\\[" + lA + "*(" + mg + ")(?:" + lA + // Operator (capture 2)
        "*([*^$|!~]?=)" + lA + // "Attribute values must be CSS identifiers [capture 5]
        // or strings [capture 3 or capture 4]"
        `*(?:'((?:\\\\.|[^\\\\'])*)'|"((?:\\\\.|[^\\\\"])*)"|(` + mg + "))|)" + lA + "*\\]", CB = ":(" + mg + `)(?:\\((('((?:\\\\.|[^\\\\'])*)'|"((?:\\\\.|[^\\\\"])*)")|((?:\\\\.|[^\\\\()[\\]]|` + zB + ")*)|.*)\\)|)", yC = new RegExp(lA + "+", "g"), LI = new RegExp("^" + lA + "+|((?:^|[^\\\\])(?:\\\\.)*)" + lA + "+$", "g"), rC = new RegExp("^" + lA + "*," + lA + "*"), $B = new RegExp("^" + lA + "*([>+~]|" + lA + ")" + lA + "*"), NC = new RegExp(lA + "|>"), kC = new RegExp(CB), nC = new RegExp("^" + mg + "$"), qI = {
          ID: new RegExp("^#(" + mg + ")"),
          CLASS: new RegExp("^\\.(" + mg + ")"),
          TAG: new RegExp("^(" + mg + "|[*])"),
          ATTR: new RegExp("^" + zB),
          PSEUDO: new RegExp("^" + CB),
          CHILD: new RegExp("^:(only|first|last|nth|nth-last)-(child|of-type)(?:\\(" + lA + "*(even|odd|(([+-]|)(\\d*)n|)" + lA + "*(?:([+-]|)" + lA + "*(\\d+)|))" + lA + "*\\)|)", "i"),
          bool: new RegExp("^(?:" + QB + ")$", "i"),
          // For use in libraries implementing .is()
          // We use this for POS matching in `select`
          needsContext: new RegExp("^" + lA + "*[>+~]|:(even|odd|eq|gt|lt|nth|first|last)(?:\\(" + lA + "*((?:-\\d)?\\d*)" + lA + "*\\)|)(?=[^-]|$)", "i")
        }, FC = /HTML$/i, RC = /^(?:input|select|textarea|button)$/i, KC = /^h\d$/i, wI = /^[^{]+\{\s*\[native \w/, YC = /^(?:#([\w-]+)|(\w+)|\.([\w-]+))$/, EB = /[+~]/, Yg = new RegExp("\\\\[\\da-fA-F]{1,6}" + lA + "?|\\\\([^\\r\\n\\f])", "g"), lg = function(y, F) {
          var J = "0x" + y.slice(1) - 65536;
          return F || // Replace a hexadecimal escape sequence with the encoded Unicode code point
          // Support: IE <=11+
          // For values outside the Basic Multilingual Plane (BMP), manually construct a
          // surrogate pair
          (J < 0 ? String.fromCharCode(J + 65536) : String.fromCharCode(J >> 10 | 55296, J & 1023 | 56320));
        }, AQ = /([\0-\x1f\x7f]|^-?\d)|^-$|[^\0-\x1f\x7f-\uFFFF\w-]/g, gQ = function(y, F) {
          return F ? y === "\0" ? "�" : y.slice(0, -1) + "\\" + y.charCodeAt(y.length - 1).toString(16) + " " : "\\" + y;
        }, IQ = function() {
          H();
        }, lC = pI(
          function(y) {
            return y.disabled === !0 && y.nodeName.toLowerCase() === "fieldset";
          },
          { dir: "parentNode", next: "legend" }
        );
        try {
          Sg.apply(
            ig = XB.call(eA.childNodes),
            eA.childNodes
          ), ig[eA.childNodes.length].nodeType;
        } catch {
          Sg = {
            apply: ig.length ? (
              // Leverage slice if possible
              function(F, J) {
                hg.apply(F, XB.call(J));
              }
            ) : (
              // Support: IE<9
              // Otherwise append directly
              function(F, J) {
                for (var m = F.length, S = 0; F[m++] = J[S++]; )
                  ;
                F.length = m - 1;
              }
            )
          };
        }
        function LA(y, F, J, m) {
          var S, b, j, P, X, aA, EA, tA = F && F.ownerDocument, NA = F ? F.nodeType : 9;
          if (J = J || [], typeof y != "string" || !y || NA !== 1 && NA !== 9 && NA !== 11)
            return J;
          if (!m && (H(F), F = F || u, cA)) {
            if (NA !== 11 && (X = YC.exec(y)))
              if (S = X[1]) {
                if (NA === 9)
                  if (j = F.getElementById(S)) {
                    if (j.id === S)
                      return J.push(j), J;
                  } else
                    return J;
                else if (tA && (j = tA.getElementById(S)) && wg(F, j) && j.id === S)
                  return J.push(j), J;
              } else {
                if (X[2])
                  return Sg.apply(J, F.getElementsByTagName(y)), J;
                if ((S = X[3]) && D.getElementsByClassName && F.getElementsByClassName)
                  return Sg.apply(J, F.getElementsByClassName(S)), J;
              }
            if (D.qsa && !tg[y + " "] && (!QA || !QA.test(y)) && // Support: IE 8 only
            // Exclude object elements
            (NA !== 1 || F.nodeName.toLowerCase() !== "object")) {
              if (EA = y, tA = F, NA === 1 && (NC.test(y) || $B.test(y))) {
                for (tA = EB.test(y) && oB(F.parentNode) || F, (tA !== F || !D.scope) && ((P = F.getAttribute("id")) ? P = P.replace(AQ, gQ) : F.setAttribute("id", P = UA)), aA = G(y), b = aA.length; b--; )
                  aA[b] = (P ? "#" + P : ":scope") + " " + HI(aA[b]);
                EA = aA.join(",");
              }
              try {
                return Sg.apply(
                  J,
                  tA.querySelectorAll(EA)
                ), J;
              } catch {
                tg(y, !0);
              } finally {
                P === UA && F.removeAttribute("id");
              }
            }
          }
          return n(y.replace(LI, "$1"), F, J, m);
        }
        function dI() {
          var y = [];
          function F(J, m) {
            return y.push(J + " ") > w.cacheLength && delete F[y.shift()], F[J + " "] = m;
          }
          return F;
        }
        function rg(y) {
          return y[UA] = !0, y;
        }
        function Ng(y) {
          var F = u.createElement("fieldset");
          try {
            return !!y(F);
          } catch {
            return !1;
          } finally {
            F.parentNode && F.parentNode.removeChild(F), F = null;
          }
        }
        function iB(y, F) {
          for (var J = y.split("|"), m = J.length; m--; )
            w.attrHandle[J[m]] = F;
        }
        function BQ(y, F) {
          var J = F && y, m = J && y.nodeType === 1 && F.nodeType === 1 && y.sourceIndex - F.sourceIndex;
          if (m)
            return m;
          if (J) {
            for (; J = J.nextSibling; )
              if (J === F)
                return -1;
          }
          return y ? 1 : -1;
        }
        function JC(y) {
          return function(F) {
            var J = F.nodeName.toLowerCase();
            return J === "input" && F.type === y;
          };
        }
        function UC(y) {
          return function(F) {
            var J = F.nodeName.toLowerCase();
            return (J === "input" || J === "button") && F.type === y;
          };
        }
        function QQ(y) {
          return function(F) {
            return "form" in F ? F.parentNode && F.disabled === !1 ? "label" in F ? "label" in F.parentNode ? F.parentNode.disabled === y : F.disabled === y : F.isDisabled === y || // Where there is no isDisabled, check manually
            /* jshint -W018 */
            F.isDisabled !== !y && lC(F) === y : F.disabled === y : "label" in F ? F.disabled === y : !1;
          };
        }
        function Og(y) {
          return rg(function(F) {
            return F = +F, rg(function(J, m) {
              for (var S, b = y([], J.length, F), j = b.length; j--; )
                J[S = b[j]] && (J[S] = !(m[S] = J[S]));
            });
          });
        }
        function oB(y) {
          return y && typeof y.getElementsByTagName < "u" && y;
        }
        D = LA.support = {}, c = LA.isXML = function(y) {
          var F = y && y.namespaceURI, J = y && (y.ownerDocument || y).documentElement;
          return !FC.test(F || J && J.nodeName || "HTML");
        }, H = LA.setDocument = function(y) {
          var F, J, m = y ? y.ownerDocument || y : eA;
          return m == u || m.nodeType !== 9 || !m.documentElement || (u = m, BA = u.documentElement, cA = !c(u), eA != u && (J = u.defaultView) && J.top !== J && (J.addEventListener ? J.addEventListener("unload", IQ, !1) : J.attachEvent && J.attachEvent("onunload", IQ)), D.scope = Ng(function(S) {
            return BA.appendChild(S).appendChild(u.createElement("div")), typeof S.querySelectorAll < "u" && !S.querySelectorAll(":scope fieldset div").length;
          }), D.attributes = Ng(function(S) {
            return S.className = "i", !S.getAttribute("className");
          }), D.getElementsByTagName = Ng(function(S) {
            return S.appendChild(u.createComment("")), !S.getElementsByTagName("*").length;
          }), D.getElementsByClassName = wI.test(u.getElementsByClassName), D.getById = Ng(function(S) {
            return BA.appendChild(S).id = UA, !u.getElementsByName || !u.getElementsByName(UA).length;
          }), D.getById ? (w.filter.ID = function(S) {
            var b = S.replace(Yg, lg);
            return function(j) {
              return j.getAttribute("id") === b;
            };
          }, w.find.ID = function(S, b) {
            if (typeof b.getElementById < "u" && cA) {
              var j = b.getElementById(S);
              return j ? [j] : [];
            }
          }) : (w.filter.ID = function(S) {
            var b = S.replace(Yg, lg);
            return function(j) {
              var P = typeof j.getAttributeNode < "u" && j.getAttributeNode("id");
              return P && P.value === b;
            };
          }, w.find.ID = function(S, b) {
            if (typeof b.getElementById < "u" && cA) {
              var j, P, X, aA = b.getElementById(S);
              if (aA) {
                if (j = aA.getAttributeNode("id"), j && j.value === S)
                  return [aA];
                for (X = b.getElementsByName(S), P = 0; aA = X[P++]; )
                  if (j = aA.getAttributeNode("id"), j && j.value === S)
                    return [aA];
              }
              return [];
            }
          }), w.find.TAG = D.getElementsByTagName ? function(S, b) {
            if (typeof b.getElementsByTagName < "u")
              return b.getElementsByTagName(S);
            if (D.qsa)
              return b.querySelectorAll(S);
          } : function(S, b) {
            var j, P = [], X = 0, aA = b.getElementsByTagName(S);
            if (S === "*") {
              for (; j = aA[X++]; )
                j.nodeType === 1 && P.push(j);
              return P;
            }
            return aA;
          }, w.find.CLASS = D.getElementsByClassName && function(S, b) {
            if (typeof b.getElementsByClassName < "u" && cA)
              return b.getElementsByClassName(S);
          }, TA = [], QA = [], (D.qsa = wI.test(u.querySelectorAll)) && (Ng(function(S) {
            var b;
            BA.appendChild(S).innerHTML = "<a id='" + UA + "'></a><select id='" + UA + "-\r\\' msallowcapture=''><option selected=''></option></select>", S.querySelectorAll("[msallowcapture^='']").length && QA.push("[*^$]=" + lA + `*(?:''|"")`), S.querySelectorAll("[selected]").length || QA.push("\\[" + lA + "*(?:value|" + QB + ")"), S.querySelectorAll("[id~=" + UA + "-]").length || QA.push("~="), b = u.createElement("input"), b.setAttribute("name", ""), S.appendChild(b), S.querySelectorAll("[name='']").length || QA.push("\\[" + lA + "*name" + lA + "*=" + lA + `*(?:''|"")`), S.querySelectorAll(":checked").length || QA.push(":checked"), S.querySelectorAll("a#" + UA + "+*").length || QA.push(".#.+[+~]"), S.querySelectorAll("\\\f"), QA.push("[\\r\\n\\f]");
          }), Ng(function(S) {
            S.innerHTML = "<a href='' disabled='disabled'></a><select disabled='disabled'><option/></select>";
            var b = u.createElement("input");
            b.setAttribute("type", "hidden"), S.appendChild(b).setAttribute("name", "D"), S.querySelectorAll("[name=d]").length && QA.push("name" + lA + "*[*^$|!~]?="), S.querySelectorAll(":enabled").length !== 2 && QA.push(":enabled", ":disabled"), BA.appendChild(S).disabled = !0, S.querySelectorAll(":disabled").length !== 2 && QA.push(":enabled", ":disabled"), S.querySelectorAll("*,:x"), QA.push(",.*:");
          })), (D.matchesSelector = wI.test(bA = BA.matches || BA.webkitMatchesSelector || BA.mozMatchesSelector || BA.oMatchesSelector || BA.msMatchesSelector)) && Ng(function(S) {
            D.disconnectedMatch = bA.call(S, "*"), bA.call(S, "[s!='']:x"), TA.push("!=", CB);
          }), QA = QA.length && new RegExp(QA.join("|")), TA = TA.length && new RegExp(TA.join("|")), F = wI.test(BA.compareDocumentPosition), wg = F || wI.test(BA.contains) ? function(S, b) {
            var j = S.nodeType === 9 ? S.documentElement : S, P = b && b.parentNode;
            return S === P || !!(P && P.nodeType === 1 && (j.contains ? j.contains(P) : S.compareDocumentPosition && S.compareDocumentPosition(P) & 16));
          } : function(S, b) {
            if (b) {
              for (; b = b.parentNode; )
                if (b === S)
                  return !0;
            }
            return !1;
          }, fg = F ? function(S, b) {
            if (S === b)
              return V = !0, 0;
            var j = !S.compareDocumentPosition - !b.compareDocumentPosition;
            return j || (j = (S.ownerDocument || S) == (b.ownerDocument || b) ? S.compareDocumentPosition(b) : (
              // Otherwise we know they are disconnected
              1
            ), j & 1 || !D.sortDetached && b.compareDocumentPosition(S) === j ? S == u || S.ownerDocument == eA && wg(eA, S) ? -1 : b == u || b.ownerDocument == eA && wg(eA, b) ? 1 : p ? xg(p, S) - xg(p, b) : 0 : j & 4 ? -1 : 1);
          } : function(S, b) {
            if (S === b)
              return V = !0, 0;
            var j, P = 0, X = S.parentNode, aA = b.parentNode, EA = [S], tA = [b];
            if (!X || !aA)
              return S == u ? -1 : b == u ? 1 : (
                /* eslint-enable eqeqeq */
                X ? -1 : aA ? 1 : p ? xg(p, S) - xg(p, b) : 0
              );
            if (X === aA)
              return BQ(S, b);
            for (j = S; j = j.parentNode; )
              EA.unshift(j);
            for (j = b; j = j.parentNode; )
              tA.unshift(j);
            for (; EA[P] === tA[P]; )
              P++;
            return P ? (
              // Do a sibling check if the nodes have a common ancestor
              BQ(EA[P], tA[P])
            ) : (
              // Otherwise nodes in our document sort first
              // Support: IE 11+, Edge 17 - 18+
              // IE/Edge sometimes throw a "Permission denied" error when strict-comparing
              // two documents; shallow comparisons work.
              /* eslint-disable eqeqeq */
              EA[P] == eA ? -1 : tA[P] == eA ? 1 : (
                /* eslint-enable eqeqeq */
                0
              )
            );
          }), u;
        }, LA.matches = function(y, F) {
          return LA(y, null, null, F);
        }, LA.matchesSelector = function(y, F) {
          if (H(y), D.matchesSelector && cA && !tg[F + " "] && (!TA || !TA.test(F)) && (!QA || !QA.test(F)))
            try {
              var J = bA.call(y, F);
              if (J || D.disconnectedMatch || // As well, disconnected nodes are said to be in a document
              // fragment in IE 9
              y.document && y.document.nodeType !== 11)
                return J;
            } catch {
              tg(F, !0);
            }
          return LA(F, u, null, [y]).length > 0;
        }, LA.contains = function(y, F) {
          return (y.ownerDocument || y) != u && H(y), wg(y, F);
        }, LA.attr = function(y, F) {
          (y.ownerDocument || y) != u && H(y);
          var J = w.attrHandle[F.toLowerCase()], m = J && ug.call(w.attrHandle, F.toLowerCase()) ? J(y, F, !cA) : void 0;
          return m !== void 0 ? m : D.attributes || !cA ? y.getAttribute(F) : (m = y.getAttributeNode(F)) && m.specified ? m.value : null;
        }, LA.escape = function(y) {
          return (y + "").replace(AQ, gQ);
        }, LA.error = function(y) {
          throw new Error("Syntax error, unrecognized expression: " + y);
        }, LA.uniqueSort = function(y) {
          var F, J = [], m = 0, S = 0;
          if (V = !D.detectDuplicates, p = !D.sortStable && y.slice(0), y.sort(fg), V) {
            for (; F = y[S++]; )
              F === y[S] && (m = J.push(S));
            for (; m--; )
              y.splice(J[m], 1);
          }
          return p = null, y;
        }, e = LA.getText = function(y) {
          var F, J = "", m = 0, S = y.nodeType;
          if (S) {
            if (S === 1 || S === 9 || S === 11) {
              if (typeof y.textContent == "string")
                return y.textContent;
              for (y = y.firstChild; y; y = y.nextSibling)
                J += e(y);
            } else if (S === 3 || S === 4)
              return y.nodeValue;
          } else
            for (; F = y[m++]; )
              J += e(F);
          return J;
        }, w = LA.selectors = {
          // Can be adjusted by the user
          cacheLength: 50,
          createPseudo: rg,
          match: qI,
          attrHandle: {},
          find: {},
          relative: {
            ">": { dir: "parentNode", first: !0 },
            " ": { dir: "parentNode" },
            "+": { dir: "previousSibling", first: !0 },
            "~": { dir: "previousSibling" }
          },
          preFilter: {
            ATTR: function(y) {
              return y[1] = y[1].replace(Yg, lg), y[3] = (y[3] || y[4] || y[5] || "").replace(Yg, lg), y[2] === "~=" && (y[3] = " " + y[3] + " "), y.slice(0, 4);
            },
            CHILD: function(y) {
              return y[1] = y[1].toLowerCase(), y[1].slice(0, 3) === "nth" ? (y[3] || LA.error(y[0]), y[4] = +(y[4] ? y[5] + (y[6] || 1) : 2 * (y[3] === "even" || y[3] === "odd")), y[5] = +(y[7] + y[8] || y[3] === "odd")) : y[3] && LA.error(y[0]), y;
            },
            PSEUDO: function(y) {
              var F, J = !y[6] && y[2];
              return qI.CHILD.test(y[0]) ? null : (y[3] ? y[2] = y[4] || y[5] || "" : J && kC.test(J) && // Get excess from tokenize (recursively)
              (F = G(J, !0)) && // advance to the next closing parenthesis
              (F = J.indexOf(")", J.length - F) - J.length) && (y[0] = y[0].slice(0, F), y[2] = J.slice(0, F)), y.slice(0, 3));
            }
          },
          filter: {
            TAG: function(y) {
              var F = y.replace(Yg, lg).toLowerCase();
              return y === "*" ? function() {
                return !0;
              } : function(J) {
                return J.nodeName && J.nodeName.toLowerCase() === F;
              };
            },
            CLASS: function(y) {
              var F = fA[y + " "];
              return F || (F = new RegExp("(^|" + lA + ")" + y + "(" + lA + "|$)")) && fA(
                y,
                function(J) {
                  return F.test(
                    typeof J.className == "string" && J.className || typeof J.getAttribute < "u" && J.getAttribute("class") || ""
                  );
                }
              );
            },
            ATTR: function(y, F, J) {
              return function(m) {
                var S = LA.attr(m, y);
                return S == null ? F === "!=" : F ? (S += "", F === "=" ? S === J : F === "!=" ? S !== J : F === "^=" ? J && S.indexOf(J) === 0 : F === "*=" ? J && S.indexOf(J) > -1 : F === "$=" ? J && S.slice(-J.length) === J : F === "~=" ? (" " + S.replace(yC, " ") + " ").indexOf(J) > -1 : F === "|=" ? S === J || S.slice(0, J.length + 1) === J + "-" : !1) : !0;
              };
            },
            CHILD: function(y, F, J, m, S) {
              var b = y.slice(0, 3) !== "nth", j = y.slice(-4) !== "last", P = F === "of-type";
              return m === 1 && S === 0 ? (
                // Shortcut for :nth-*(n)
                function(X) {
                  return !!X.parentNode;
                }
              ) : function(X, aA, EA) {
                var tA, NA, qA, wA, WA, Ag, sg = b !== j ? "nextSibling" : "previousSibling", pA = X.parentNode, tI = P && X.nodeName.toLowerCase(), sI = !EA && !P, eg = !1;
                if (pA) {
                  if (b) {
                    for (; sg; ) {
                      for (wA = X; wA = wA[sg]; )
                        if (P ? wA.nodeName.toLowerCase() === tI : wA.nodeType === 1)
                          return !1;
                      Ag = sg = y === "only" && !Ag && "nextSibling";
                    }
                    return !0;
                  }
                  if (Ag = [j ? pA.firstChild : pA.lastChild], j && sI) {
                    for (wA = pA, qA = wA[UA] || (wA[UA] = {}), NA = qA[wA.uniqueID] || (qA[wA.uniqueID] = {}), tA = NA[y] || [], WA = tA[0] === Eg && tA[1], eg = WA && tA[2], wA = WA && pA.childNodes[WA]; wA = ++WA && wA && wA[sg] || // Fallback to seeking `elem` from the start
                    (eg = WA = 0) || Ag.pop(); )
                      if (wA.nodeType === 1 && ++eg && wA === X) {
                        NA[y] = [Eg, WA, eg];
                        break;
                      }
                  } else if (sI && (wA = X, qA = wA[UA] || (wA[UA] = {}), NA = qA[wA.uniqueID] || (qA[wA.uniqueID] = {}), tA = NA[y] || [], WA = tA[0] === Eg && tA[1], eg = WA), eg === !1)
                    for (; (wA = ++WA && wA && wA[sg] || (eg = WA = 0) || Ag.pop()) && !((P ? wA.nodeName.toLowerCase() === tI : wA.nodeType === 1) && ++eg && (sI && (qA = wA[UA] || (wA[UA] = {}), NA = qA[wA.uniqueID] || (qA[wA.uniqueID] = {}), NA[y] = [Eg, eg]), wA === X)); )
                      ;
                  return eg -= S, eg === m || eg % m === 0 && eg / m >= 0;
                }
              };
            },
            PSEUDO: function(y, F) {
              var J, m = w.pseudos[y] || w.setFilters[y.toLowerCase()] || LA.error("unsupported pseudo: " + y);
              return m[UA] ? m(F) : m.length > 1 ? (J = [y, y, "", F], w.setFilters.hasOwnProperty(y.toLowerCase()) ? rg(function(S, b) {
                for (var j, P = m(S, F), X = P.length; X--; )
                  j = xg(S, P[X]), S[j] = !(b[j] = P[X]);
              }) : function(S) {
                return m(S, 0, J);
              }) : m;
            }
          },
          pseudos: {
            // Potentially complex pseudos
            not: rg(function(y) {
              var F = [], J = [], m = R(y.replace(LI, "$1"));
              return m[UA] ? rg(function(S, b, j, P) {
                for (var X, aA = m(S, null, P, []), EA = S.length; EA--; )
                  (X = aA[EA]) && (S[EA] = !(b[EA] = X));
              }) : function(S, b, j) {
                return F[0] = S, m(F, null, j, J), F[0] = null, !J.pop();
              };
            }),
            has: rg(function(y) {
              return function(F) {
                return LA(y, F).length > 0;
              };
            }),
            contains: rg(function(y) {
              return y = y.replace(Yg, lg), function(F) {
                return (F.textContent || e(F)).indexOf(y) > -1;
              };
            }),
            // "Whether an element is represented by a :lang() selector
            // is based solely on the element's language value
            // being equal to the identifier C,
            // or beginning with the identifier C immediately followed by "-".
            // The matching of C against the element's language value is performed case-insensitively.
            // The identifier C does not have to be a valid language name."
            // http://www.w3.org/TR/selectors/#lang-pseudo
            lang: rg(function(y) {
              return nC.test(y || "") || LA.error("unsupported lang: " + y), y = y.replace(Yg, lg).toLowerCase(), function(F) {
                var J;
                do
                  if (J = cA ? F.lang : F.getAttribute("xml:lang") || F.getAttribute("lang"))
                    return J = J.toLowerCase(), J === y || J.indexOf(y + "-") === 0;
                while ((F = F.parentNode) && F.nodeType === 1);
                return !1;
              };
            }),
            // Miscellaneous
            target: function(y) {
              var F = Q.location && Q.location.hash;
              return F && F.slice(1) === y.id;
            },
            root: function(y) {
              return y === BA;
            },
            focus: function(y) {
              return y === u.activeElement && (!u.hasFocus || u.hasFocus()) && !!(y.type || y.href || ~y.tabIndex);
            },
            // Boolean properties
            enabled: QQ(!1),
            disabled: QQ(!0),
            checked: function(y) {
              var F = y.nodeName.toLowerCase();
              return F === "input" && !!y.checked || F === "option" && !!y.selected;
            },
            selected: function(y) {
              return y.parentNode && y.parentNode.selectedIndex, y.selected === !0;
            },
            // Contents
            empty: function(y) {
              for (y = y.firstChild; y; y = y.nextSibling)
                if (y.nodeType < 6)
                  return !1;
              return !0;
            },
            parent: function(y) {
              return !w.pseudos.empty(y);
            },
            // Element/input types
            header: function(y) {
              return KC.test(y.nodeName);
            },
            input: function(y) {
              return RC.test(y.nodeName);
            },
            button: function(y) {
              var F = y.nodeName.toLowerCase();
              return F === "input" && y.type === "button" || F === "button";
            },
            text: function(y) {
              var F;
              return y.nodeName.toLowerCase() === "input" && y.type === "text" && // Support: IE<8
              // New HTML5 attribute values (e.g., "search") appear with elem.type === "text"
              ((F = y.getAttribute("type")) == null || F.toLowerCase() === "text");
            },
            // Position-in-collection
            first: Og(function() {
              return [0];
            }),
            last: Og(function(y, F) {
              return [F - 1];
            }),
            eq: Og(function(y, F, J) {
              return [J < 0 ? J + F : J];
            }),
            even: Og(function(y, F) {
              for (var J = 0; J < F; J += 2)
                y.push(J);
              return y;
            }),
            odd: Og(function(y, F) {
              for (var J = 1; J < F; J += 2)
                y.push(J);
              return y;
            }),
            lt: Og(function(y, F, J) {
              for (var m = J < 0 ? J + F : J > F ? F : J; --m >= 0; )
                y.push(m);
              return y;
            }),
            gt: Og(function(y, F, J) {
              for (var m = J < 0 ? J + F : J; ++m < F; )
                y.push(m);
              return y;
            })
          }
        }, w.pseudos.nth = w.pseudos.eq;
        for (i in { radio: !0, checkbox: !0, file: !0, password: !0, image: !0 })
          w.pseudos[i] = JC(i);
        for (i in { submit: !0, reset: !0 })
          w.pseudos[i] = UC(i);
        function CQ() {
        }
        CQ.prototype = w.filters = w.pseudos, w.setFilters = new CQ(), G = LA.tokenize = function(y, F) {
          var J, m, S, b, j, P, X, aA = aI[y + " "];
          if (aA)
            return F ? 0 : aA.slice(0);
          for (j = y, P = [], X = w.preFilter; j; ) {
            (!J || (m = rC.exec(j))) && (m && (j = j.slice(m[0].length) || j), P.push(S = [])), J = !1, (m = $B.exec(j)) && (J = m.shift(), S.push({
              value: J,
              // Cast descendant combinators to space
              type: m[0].replace(LI, " ")
            }), j = j.slice(J.length));
            for (b in w.filter)
              (m = qI[b].exec(j)) && (!X[b] || (m = X[b](m))) && (J = m.shift(), S.push({
                value: J,
                type: b,
                matches: m
              }), j = j.slice(J.length));
            if (!J)
              break;
          }
          return F ? j.length : j ? LA.error(y) : (
            // Cache the tokens
            aI(y, P).slice(0)
          );
        };
        function HI(y) {
          for (var F = 0, J = y.length, m = ""; F < J; F++)
            m += y[F].value;
          return m;
        }
        function pI(y, F, J) {
          var m = F.dir, S = F.next, b = S || m, j = J && b === "parentNode", P = YA++;
          return F.first ? (
            // Check against closest ancestor/preceding element
            function(X, aA, EA) {
              for (; X = X[m]; )
                if (X.nodeType === 1 || j)
                  return y(X, aA, EA);
              return !1;
            }
          ) : (
            // Check against all ancestor/preceding elements
            function(X, aA, EA) {
              var tA, NA, qA, wA = [Eg, P];
              if (EA) {
                for (; X = X[m]; )
                  if ((X.nodeType === 1 || j) && y(X, aA, EA))
                    return !0;
              } else
                for (; X = X[m]; )
                  if (X.nodeType === 1 || j)
                    if (qA = X[UA] || (X[UA] = {}), NA = qA[X.uniqueID] || (qA[X.uniqueID] = {}), S && S === X.nodeName.toLowerCase())
                      X = X[m] || X;
                    else {
                      if ((tA = NA[b]) && tA[0] === Eg && tA[1] === P)
                        return wA[2] = tA[2];
                      if (NA[b] = wA, wA[2] = y(X, aA, EA))
                        return !0;
                    }
              return !1;
            }
          );
        }
        function DB(y) {
          return y.length > 1 ? function(F, J, m) {
            for (var S = y.length; S--; )
              if (!y[S](F, J, m))
                return !1;
            return !0;
          } : y[0];
        }
        function SC(y, F, J) {
          for (var m = 0, S = F.length; m < S; m++)
            LA(y, F[m], J);
          return J;
        }
        function fI(y, F, J, m, S) {
          for (var b, j = [], P = 0, X = y.length, aA = F != null; P < X; P++)
            (b = y[P]) && (!J || J(b, m, S)) && (j.push(b), aA && F.push(P));
          return j;
        }
        function aB(y, F, J, m, S, b) {
          return m && !m[UA] && (m = aB(m)), S && !S[UA] && (S = aB(S, b)), rg(function(j, P, X, aA) {
            var EA, tA, NA, qA = [], wA = [], WA = P.length, Ag = j || SC(
              F || "*",
              X.nodeType ? [X] : X,
              []
            ), sg = y && (j || !F) ? fI(Ag, qA, y, X, aA) : Ag, pA = J ? (
              // If we have a postFinder, or filtered seed, or non-seed postFilter or preexisting results,
              S || (j ? y : WA || m) ? (
                // ...intermediate processing is necessary
                []
              ) : (
                // ...otherwise use results directly
                P
              )
            ) : sg;
            if (J && J(sg, pA, X, aA), m)
              for (EA = fI(pA, wA), m(EA, [], X, aA), tA = EA.length; tA--; )
                (NA = EA[tA]) && (pA[wA[tA]] = !(sg[wA[tA]] = NA));
            if (j) {
              if (S || y) {
                if (S) {
                  for (EA = [], tA = pA.length; tA--; )
                    (NA = pA[tA]) && EA.push(sg[tA] = NA);
                  S(null, pA = [], EA, aA);
                }
                for (tA = pA.length; tA--; )
                  (NA = pA[tA]) && (EA = S ? xg(j, NA) : qA[tA]) > -1 && (j[EA] = !(P[EA] = NA));
              }
            } else
              pA = fI(
                pA === P ? pA.splice(WA, pA.length) : pA
              ), S ? S(null, P, pA, aA) : Sg.apply(P, pA);
          });
        }
        function wB(y) {
          for (var F, J, m, S = y.length, b = w.relative[y[0].type], j = b || w.relative[" "], P = b ? 1 : 0, X = pI(function(tA) {
            return tA === F;
          }, j, !0), aA = pI(function(tA) {
            return xg(F, tA) > -1;
          }, j, !0), EA = [function(tA, NA, qA) {
            var wA = !b && (qA || NA !== L) || ((F = NA).nodeType ? X(tA, NA, qA) : aA(tA, NA, qA));
            return F = null, wA;
          }]; P < S; P++)
            if (J = w.relative[y[P].type])
              EA = [pI(DB(EA), J)];
            else {
              if (J = w.filter[y[P].type].apply(null, y[P].matches), J[UA]) {
                for (m = ++P; m < S && !w.relative[y[m].type]; m++)
                  ;
                return aB(
                  P > 1 && DB(EA),
                  P > 1 && HI(
                    // If the preceding token was a descendant combinator, insert an implicit any-element `*`
                    y.slice(0, P - 1).concat({ value: y[P - 2].type === " " ? "*" : "" })
                  ).replace(LI, "$1"),
                  J,
                  P < m && wB(y.slice(P, m)),
                  m < S && wB(y = y.slice(m)),
                  m < S && HI(y)
                );
              }
              EA.push(J);
            }
          return DB(EA);
        }
        function LC(y, F) {
          var J = F.length > 0, m = y.length > 0, S = function(b, j, P, X, aA) {
            var EA, tA, NA, qA = 0, wA = "0", WA = b && [], Ag = [], sg = L, pA = b || m && w.find.TAG("*", aA), tI = Eg += sg == null ? 1 : Math.random() || 0.1, sI = pA.length;
            for (aA && (L = j == u || j || aA); wA !== sI && (EA = pA[wA]) != null; wA++) {
              if (m && EA) {
                for (tA = 0, !j && EA.ownerDocument != u && (H(EA), P = !cA); NA = y[tA++]; )
                  if (NA(EA, j || u, P)) {
                    X.push(EA);
                    break;
                  }
                aA && (Eg = tI);
              }
              J && ((EA = !NA && EA) && qA--, b && WA.push(EA));
            }
            if (qA += wA, J && wA !== qA) {
              for (tA = 0; NA = F[tA++]; )
                NA(WA, Ag, j, P);
              if (b) {
                if (qA > 0)
                  for (; wA--; )
                    WA[wA] || Ag[wA] || (Ag[wA] = Ug.call(X));
                Ag = fI(Ag);
              }
              Sg.apply(X, Ag), aA && !b && Ag.length > 0 && qA + F.length > 1 && LA.uniqueSort(X);
            }
            return aA && (Eg = tI, L = sg), WA;
          };
          return J ? rg(S) : S;
        }
        return R = LA.compile = function(y, F) {
          var J, m = [], S = [], b = SI[y + " "];
          if (!b) {
            for (F || (F = G(y)), J = F.length; J--; )
              b = wB(F[J]), b[UA] ? m.push(b) : S.push(b);
            b = SI(
              y,
              LC(S, m)
            ), b.selector = y;
          }
          return b;
        }, n = LA.select = function(y, F, J, m) {
          var S, b, j, P, X, aA = typeof y == "function" && y, EA = !m && G(y = aA.selector || y);
          if (J = J || [], EA.length === 1) {
            if (b = EA[0] = EA[0].slice(0), b.length > 2 && (j = b[0]).type === "ID" && F.nodeType === 9 && cA && w.relative[b[1].type]) {
              if (F = (w.find.ID(j.matches[0].replace(Yg, lg), F) || [])[0], F)
                aA && (F = F.parentNode);
              else
                return J;
              y = y.slice(b.shift().value.length);
            }
            for (S = qI.needsContext.test(y) ? 0 : b.length; S-- && (j = b[S], !w.relative[P = j.type]); )
              if ((X = w.find[P]) && (m = X(
                j.matches[0].replace(Yg, lg),
                EB.test(b[0].type) && oB(F.parentNode) || F
              ))) {
                if (b.splice(S, 1), y = m.length && HI(b), !y)
                  return Sg.apply(J, m), J;
                break;
              }
          }
          return (aA || R(y, EA))(
            m,
            F,
            !cA,
            J,
            !F || EB.test(y) && oB(F.parentNode) || F
          ), J;
        }, D.sortStable = UA.split("").sort(fg).join("") === UA, D.detectDuplicates = !!V, H(), D.sortDetached = Ng(function(y) {
          return y.compareDocumentPosition(u.createElement("fieldset")) & 1;
        }), Ng(function(y) {
          return y.innerHTML = "<a href='#'></a>", y.firstChild.getAttribute("href") === "#";
        }) || iB("type|href|height|width", function(y, F, J) {
          if (!J)
            return y.getAttribute(F, F.toLowerCase() === "type" ? 1 : 2);
        }), (!D.attributes || !Ng(function(y) {
          return y.innerHTML = "<input/>", y.firstChild.setAttribute("value", ""), y.firstChild.getAttribute("value") === "";
        })) && iB("value", function(y, F, J) {
          if (!J && y.nodeName.toLowerCase() === "input")
            return y.defaultValue;
        }), Ng(function(y) {
          return y.getAttribute("disabled") == null;
        }) || iB(QB, function(y, F, J) {
          var m;
          if (!J)
            return y[F] === !0 ? F.toLowerCase() : (m = y.getAttributeNode(F)) && m.specified ? m.value : null;
        }), LA;
      }(A)
    );
    a.find = IA, a.expr = IA.selectors, a.expr[":"] = a.expr.pseudos, a.uniqueSort = a.unique = IA.uniqueSort, a.text = IA.getText, a.isXMLDoc = IA.isXML, a.contains = IA.contains, a.escapeSelector = IA.escape;
    var z = function(Q, i, D) {
      for (var w = [], e = D !== void 0; (Q = Q[i]) && Q.nodeType !== 9; )
        if (Q.nodeType === 1) {
          if (e && a(Q).is(D))
            break;
          w.push(Q);
        }
      return w;
    }, oA = function(Q, i) {
      for (var D = []; Q; Q = Q.nextSibling)
        Q.nodeType === 1 && Q !== i && D.push(Q);
      return D;
    }, hA = a.expr.match.needsContext;
    function yA(Q, i) {
      return Q.nodeName && Q.nodeName.toLowerCase() === i.toLowerCase();
    }
    var rA = /^<([a-z][^\/\0>:\x20\t\r\n\f]*)[\x20\t\r\n\f]*\/?>(?:<\/\1>|)$/i;
    function gA(Q, i, D) {
      return l(i) ? a.grep(Q, function(w, e) {
        return !!i.call(w, e, w) !== D;
      }) : i.nodeType ? a.grep(Q, function(w) {
        return w === i !== D;
      }) : typeof i != "string" ? a.grep(Q, function(w) {
        return s.call(i, w) > -1 !== D;
      }) : a.filter(i, Q, D);
    }
    a.filter = function(Q, i, D) {
      var w = i[0];
      return D && (Q = ":not(" + Q + ")"), i.length === 1 && w.nodeType === 1 ? a.find.matchesSelector(w, Q) ? [w] : [] : a.find.matches(Q, a.grep(i, function(e) {
        return e.nodeType === 1;
      }));
    }, a.fn.extend({
      find: function(Q) {
        var i, D, w = this.length, e = this;
        if (typeof Q != "string")
          return this.pushStack(a(Q).filter(function() {
            for (i = 0; i < w; i++)
              if (a.contains(e[i], this))
                return !0;
          }));
        for (D = this.pushStack([]), i = 0; i < w; i++)
          a.find(Q, e[i], D);
        return w > 1 ? a.uniqueSort(D) : D;
      },
      filter: function(Q) {
        return this.pushStack(gA(this, Q || [], !1));
      },
      not: function(Q) {
        return this.pushStack(gA(this, Q || [], !0));
      },
      is: function(Q) {
        return !!gA(
          this,
          // If this is a positional/relative selector, check membership in the returned set
          // so $("p:first").is("p:last") won't return true for a doc with two "p".
          typeof Q == "string" && hA.test(Q) ? a(Q) : Q || [],
          !1
        ).length;
      }
    });
    var nA, og = /^(?:\s*(<[\w\W]+>)[^>]*|#([\w-]+))$/, PA = a.fn.init = function(Q, i, D) {
      var w, e;
      if (!Q)
        return this;
      if (D = D || nA, typeof Q == "string")
        if (Q[0] === "<" && Q[Q.length - 1] === ">" && Q.length >= 3 ? w = [null, Q, null] : w = og.exec(Q), w && (w[1] || !i))
          if (w[1]) {
            if (i = i instanceof a ? i[0] : i, a.merge(this, a.parseHTML(
              w[1],
              i && i.nodeType ? i.ownerDocument || i : q,
              !0
            )), rA.test(w[1]) && a.isPlainObject(i))
              for (w in i)
                l(this[w]) ? this[w](i[w]) : this.attr(w, i[w]);
            return this;
          } else
            return e = q.getElementById(w[2]), e && (this[0] = e, this.length = 1), this;
        else
          return !i || i.jquery ? (i || D).find(Q) : this.constructor(i).find(Q);
      else {
        if (Q.nodeType)
          return this[0] = Q, this.length = 1, this;
        if (l(Q))
          return D.ready !== void 0 ? D.ready(Q) : (
            // Execute immediately if ready is not present
            Q(a)
          );
      }
      return a.makeArray(Q, this);
    };
    PA.prototype = a.fn, nA = a(q);
    var gg = /^(?:parents|prev(?:Until|All))/, XA = {
      children: !0,
      contents: !0,
      next: !0,
      prev: !0
    };
    a.fn.extend({
      has: function(Q) {
        var i = a(Q, this), D = i.length;
        return this.filter(function() {
          for (var w = 0; w < D; w++)
            if (a.contains(this, i[w]))
              return !0;
        });
      },
      closest: function(Q, i) {
        var D, w = 0, e = this.length, c = [], G = typeof Q != "string" && a(Q);
        if (!hA.test(Q)) {
          for (; w < e; w++)
            for (D = this[w]; D && D !== i; D = D.parentNode)
              if (D.nodeType < 11 && (G ? G.index(D) > -1 : (
                // Don't pass non-elements to Sizzle
                D.nodeType === 1 && a.find.matchesSelector(D, Q)
              ))) {
                c.push(D);
                break;
              }
        }
        return this.pushStack(c.length > 1 ? a.uniqueSort(c) : c);
      },
      // Determine the position of an element within the set
      index: function(Q) {
        return Q ? typeof Q == "string" ? s.call(a(Q), this[0]) : s.call(
          this,
          // If it receives a jQuery object, the first element is used
          Q.jquery ? Q[0] : Q
        ) : this[0] && this[0].parentNode ? this.first().prevAll().length : -1;
      },
      add: function(Q, i) {
        return this.pushStack(
          a.uniqueSort(
            a.merge(this.get(), a(Q, i))
          )
        );
      },
      addBack: function(Q) {
        return this.add(
          Q == null ? this.prevObject : this.prevObject.filter(Q)
        );
      }
    });
    function OA(Q, i) {
      for (; (Q = Q[i]) && Q.nodeType !== 1; )
        ;
      return Q;
    }
    a.each({
      parent: function(Q) {
        var i = Q.parentNode;
        return i && i.nodeType !== 11 ? i : null;
      },
      parents: function(Q) {
        return z(Q, "parentNode");
      },
      parentsUntil: function(Q, i, D) {
        return z(Q, "parentNode", D);
      },
      next: function(Q) {
        return OA(Q, "nextSibling");
      },
      prev: function(Q) {
        return OA(Q, "previousSibling");
      },
      nextAll: function(Q) {
        return z(Q, "nextSibling");
      },
      prevAll: function(Q) {
        return z(Q, "previousSibling");
      },
      nextUntil: function(Q, i, D) {
        return z(Q, "nextSibling", D);
      },
      prevUntil: function(Q, i, D) {
        return z(Q, "previousSibling", D);
      },
      siblings: function(Q) {
        return oA((Q.parentNode || {}).firstChild, Q);
      },
      children: function(Q) {
        return oA(Q.firstChild);
      },
      contents: function(Q) {
        return Q.contentDocument != null && // Support: IE 11+
        // <object> elements with no `data` attribute has an object
        // `contentDocument` with a `null` prototype.
        I(Q.contentDocument) ? Q.contentDocument : (yA(Q, "template") && (Q = Q.content || Q), a.merge([], Q.childNodes));
      }
    }, function(Q, i) {
      a.fn[Q] = function(D, w) {
        var e = a.map(this, i, D);
        return Q.slice(-5) !== "Until" && (w = D), w && typeof w == "string" && (e = a.filter(w, e)), this.length > 1 && (XA[Q] || a.uniqueSort(e), gg.test(Q) && e.reverse()), this.pushStack(e);
      };
    });
    var Dg = /[^\x20\t\r\n\f]+/g;
    function RI(Q) {
      var i = {};
      return a.each(Q.match(Dg) || [], function(D, w) {
        i[w] = !0;
      }), i;
    }
    a.Callbacks = function(Q) {
      Q = typeof Q == "string" ? RI(Q) : a.extend({}, Q);
      var i, D, w, e, c = [], G = [], R = -1, n = function() {
        for (e = e || Q.once, w = i = !0; G.length; R = -1)
          for (D = G.shift(); ++R < c.length; )
            c[R].apply(D[0], D[1]) === !1 && Q.stopOnFalse && (R = c.length, D = !1);
        Q.memory || (D = !1), i = !1, e && (D ? c = [] : c = "");
      }, L = {
        // Add a callback or a collection of callbacks to the list
        add: function() {
          return c && (D && !i && (R = c.length - 1, G.push(D)), function p(V) {
            a.each(V, function(H, u) {
              l(u) ? (!Q.unique || !L.has(u)) && c.push(u) : u && u.length && v(u) !== "string" && p(u);
            });
          }(arguments), D && !i && n()), this;
        },
        // Remove a callback from the list
        remove: function() {
          return a.each(arguments, function(p, V) {
            for (var H; (H = a.inArray(V, c, H)) > -1; )
              c.splice(H, 1), H <= R && R--;
          }), this;
        },
        // Check if a given callback is in the list.
        // If no argument is given, return whether or not list has callbacks attached.
        has: function(p) {
          return p ? a.inArray(p, c) > -1 : c.length > 0;
        },
        // Remove all callbacks from the list
        empty: function() {
          return c && (c = []), this;
        },
        // Disable .fire and .add
        // Abort any current/pending executions
        // Clear all callbacks and values
        disable: function() {
          return e = G = [], c = D = "", this;
        },
        disabled: function() {
          return !c;
        },
        // Disable .fire
        // Also disable .add unless we have memory (since it would have no effect)
        // Abort any pending executions
        lock: function() {
          return e = G = [], !D && !i && (c = D = ""), this;
        },
        locked: function() {
          return !!e;
        },
        // Call all callbacks with the given context and arguments
        fireWith: function(p, V) {
          return e || (V = V || [], V = [p, V.slice ? V.slice() : V], G.push(V), i || n()), this;
        },
        // Call all the callbacks with the given arguments
        fire: function() {
          return L.fireWith(this, arguments), this;
        },
        // To know if the callbacks have already been called at least once
        fired: function() {
          return !!w;
        }
      };
      return L;
    };
    function Jg(Q) {
      return Q;
    }
    function sA(Q) {
      throw Q;
    }
    function KA(Q, i, D, w) {
      var e;
      try {
        Q && l(e = Q.promise) ? e.call(Q).done(i).fail(D) : Q && l(e = Q.then) ? e.call(Q, i, D) : i.apply(void 0, [Q].slice(w));
      } catch (c) {
        D.apply(void 0, [c]);
      }
    }
    a.extend({
      Deferred: function(Q) {
        var i = [
          // action, add listener, callbacks,
          // ... .then handlers, argument index, [final state]
          [
            "notify",
            "progress",
            a.Callbacks("memory"),
            a.Callbacks("memory"),
            2
          ],
          [
            "resolve",
            "done",
            a.Callbacks("once memory"),
            a.Callbacks("once memory"),
            0,
            "resolved"
          ],
          [
            "reject",
            "fail",
            a.Callbacks("once memory"),
            a.Callbacks("once memory"),
            1,
            "rejected"
          ]
        ], D = "pending", w = {
          state: function() {
            return D;
          },
          always: function() {
            return e.done(arguments).fail(arguments), this;
          },
          catch: function(c) {
            return w.then(null, c);
          },
          // Keep pipe for back-compat
          pipe: function() {
            var c = arguments;
            return a.Deferred(function(G) {
              a.each(i, function(R, n) {
                var L = l(c[n[4]]) && c[n[4]];
                e[n[1]](function() {
                  var p = L && L.apply(this, arguments);
                  p && l(p.promise) ? p.promise().progress(G.notify).done(G.resolve).fail(G.reject) : G[n[0] + "With"](
                    this,
                    L ? [p] : arguments
                  );
                });
              }), c = null;
            }).promise();
          },
          then: function(c, G, R) {
            var n = 0;
            function L(p, V, H, u) {
              return function() {
                var BA = this, cA = arguments, QA = function() {
                  var bA, wg;
                  if (!(p < n)) {
                    if (bA = H.apply(BA, cA), bA === V.promise())
                      throw new TypeError("Thenable self-resolution");
                    wg = bA && // Support: Promises/A+ section 2.3.4
                    // https://promisesaplus.com/#point-64
                    // Only check objects and functions for thenability
                    (typeof bA == "object" || typeof bA == "function") && bA.then, l(wg) ? u ? wg.call(
                      bA,
                      L(n, V, Jg, u),
                      L(n, V, sA, u)
                    ) : (n++, wg.call(
                      bA,
                      L(n, V, Jg, u),
                      L(n, V, sA, u),
                      L(
                        n,
                        V,
                        Jg,
                        V.notifyWith
                      )
                    )) : (H !== Jg && (BA = void 0, cA = [bA]), (u || V.resolveWith)(BA, cA));
                  }
                }, TA = u ? QA : function() {
                  try {
                    QA();
                  } catch (bA) {
                    a.Deferred.exceptionHook && a.Deferred.exceptionHook(
                      bA,
                      TA.stackTrace
                    ), p + 1 >= n && (H !== sA && (BA = void 0, cA = [bA]), V.rejectWith(BA, cA));
                  }
                };
                p ? TA() : (a.Deferred.getStackHook && (TA.stackTrace = a.Deferred.getStackHook()), A.setTimeout(TA));
              };
            }
            return a.Deferred(function(p) {
              i[0][3].add(
                L(
                  0,
                  p,
                  l(R) ? R : Jg,
                  p.notifyWith
                )
              ), i[1][3].add(
                L(
                  0,
                  p,
                  l(c) ? c : Jg
                )
              ), i[2][3].add(
                L(
                  0,
                  p,
                  l(G) ? G : sA
                )
              );
            }).promise();
          },
          // Get a promise for this deferred
          // If obj is provided, the promise aspect is added to the object
          promise: function(c) {
            return c != null ? a.extend(c, w) : w;
          }
        }, e = {};
        return a.each(i, function(c, G) {
          var R = G[2], n = G[5];
          w[G[1]] = R.add, n && R.add(
            function() {
              D = n;
            },
            // rejected_callbacks.disable
            // fulfilled_callbacks.disable
            i[3 - c][2].disable,
            // rejected_handlers.disable
            // fulfilled_handlers.disable
            i[3 - c][3].disable,
            // progress_callbacks.lock
            i[0][2].lock,
            // progress_handlers.lock
            i[0][3].lock
          ), R.add(G[3].fire), e[G[0]] = function() {
            return e[G[0] + "With"](this === e ? void 0 : this, arguments), this;
          }, e[G[0] + "With"] = R.fireWith;
        }), w.promise(e), Q && Q.call(e, e), e;
      },
      // Deferred helper
      when: function(Q) {
        var i = arguments.length, D = i, w = Array(D), e = E.call(arguments), c = a.Deferred(), G = function(R) {
          return function(n) {
            w[R] = this, e[R] = arguments.length > 1 ? E.call(arguments) : n, --i || c.resolveWith(w, e);
          };
        };
        if (i <= 1 && (KA(
          Q,
          c.done(G(D)).resolve,
          c.reject,
          !i
        ), c.state() === "pending" || l(e[D] && e[D].then)))
          return c.then();
        for (; D--; )
          KA(e[D], G(D), c.reject);
        return c.promise();
      }
    });
    var xA = /^(Eval|Internal|Range|Reference|Syntax|Type|URI)Error$/;
    a.Deferred.exceptionHook = function(Q, i) {
      A.console && A.console.warn && Q && xA.test(Q.name) && A.console.warn("jQuery.Deferred exception: " + Q.message, Q.stack, i);
    }, a.readyException = function(Q) {
      A.setTimeout(function() {
        throw Q;
      });
    };
    var vA = a.Deferred();
    a.fn.ready = function(Q) {
      return vA.then(Q).catch(function(i) {
        a.readyException(i);
      }), this;
    }, a.extend({
      // Is the DOM ready to be used? Set to true once it occurs.
      isReady: !1,
      // A counter to track how many items to wait for before
      // the ready event fires. See trac-6781
      readyWait: 1,
      // Handle when the DOM is ready
      ready: function(Q) {
        (Q === !0 ? --a.readyWait : a.isReady) || (a.isReady = !0, !(Q !== !0 && --a.readyWait > 0) && vA.resolveWith(q, [a]));
      }
    }), a.ready.then = vA.then;
    function ZA() {
      q.removeEventListener("DOMContentLoaded", ZA), A.removeEventListener("load", ZA), a.ready();
    }
    q.readyState === "complete" || q.readyState !== "loading" && !q.documentElement.doScroll ? A.setTimeout(a.ready) : (q.addEventListener("DOMContentLoaded", ZA), A.addEventListener("load", ZA));
    var HA = function(Q, i, D, w, e, c, G) {
      var R = 0, n = Q.length, L = D == null;
      if (v(D) === "object") {
        e = !0;
        for (R in D)
          HA(Q, i, R, D[R], !0, c, G);
      } else if (w !== void 0 && (e = !0, l(w) || (G = !0), L && (G ? (i.call(Q, w), i = null) : (L = i, i = function(p, V, H) {
        return L.call(a(p), H);
      })), i))
        for (; R < n; R++)
          i(
            Q[R],
            D,
            G ? w : w.call(Q[R], R, i(Q[R], D))
          );
      return e ? Q : L ? i.call(Q) : n ? i(Q[0], D) : c;
    }, Rg = /^-ms-/, Ig = /-([a-z])/g;
    function qg(Q, i) {
      return i.toUpperCase();
    }
    function ag(Q) {
      return Q.replace(Rg, "ms-").replace(Ig, qg);
    }
    var II = function(Q) {
      return Q.nodeType === 1 || Q.nodeType === 9 || !+Q.nodeType;
    };
    function BI() {
      this.expando = a.expando + BI.uid++;
    }
    BI.uid = 1, BI.prototype = {
      cache: function(Q) {
        var i = Q[this.expando];
        return i || (i = {}, II(Q) && (Q.nodeType ? Q[this.expando] = i : Object.defineProperty(Q, this.expando, {
          value: i,
          configurable: !0
        }))), i;
      },
      set: function(Q, i, D) {
        var w, e = this.cache(Q);
        if (typeof i == "string")
          e[ag(i)] = D;
        else
          for (w in i)
            e[ag(w)] = i[w];
        return e;
      },
      get: function(Q, i) {
        return i === void 0 ? this.cache(Q) : (
          // Always use camelCase key (gh-2257)
          Q[this.expando] && Q[this.expando][ag(i)]
        );
      },
      access: function(Q, i, D) {
        return i === void 0 || i && typeof i == "string" && D === void 0 ? this.get(Q, i) : (this.set(Q, i, D), D !== void 0 ? D : i);
      },
      remove: function(Q, i) {
        var D, w = Q[this.expando];
        if (w !== void 0) {
          if (i !== void 0)
            for (Array.isArray(i) ? i = i.map(ag) : (i = ag(i), i = i in w ? [i] : i.match(Dg) || []), D = i.length; D--; )
              delete w[i[D]];
          (i === void 0 || a.isEmptyObject(w)) && (Q.nodeType ? Q[this.expando] = void 0 : delete Q[this.expando]);
        }
      },
      hasData: function(Q) {
        var i = Q[this.expando];
        return i !== void 0 && !a.isEmptyObject(i);
      }
    };
    var DA = new BI(), Bg = new BI(), SQ = /^(?:\{[\w\W]*\}|\[[\w\W]*\])$/, LQ = /[A-Z]/g;
    function qQ(Q) {
      return Q === "true" ? !0 : Q === "false" ? !1 : Q === "null" ? null : Q === +Q + "" ? +Q : SQ.test(Q) ? JSON.parse(Q) : Q;
    }
    function yB(Q, i, D) {
      var w;
      if (D === void 0 && Q.nodeType === 1)
        if (w = "data-" + i.replace(LQ, "-$&").toLowerCase(), D = Q.getAttribute(w), typeof D == "string") {
          try {
            D = qQ(D);
          } catch {
          }
          Bg.set(Q, i, D);
        } else
          D = void 0;
      return D;
    }
    a.extend({
      hasData: function(Q) {
        return Bg.hasData(Q) || DA.hasData(Q);
      },
      data: function(Q, i, D) {
        return Bg.access(Q, i, D);
      },
      removeData: function(Q, i) {
        Bg.remove(Q, i);
      },
      // TODO: Now that all calls to _data and _removeData have been replaced
      // with direct calls to dataPriv methods, these can be deprecated.
      _data: function(Q, i, D) {
        return DA.access(Q, i, D);
      },
      _removeData: function(Q, i) {
        DA.remove(Q, i);
      }
    }), a.fn.extend({
      data: function(Q, i) {
        var D, w, e, c = this[0], G = c && c.attributes;
        if (Q === void 0) {
          if (this.length && (e = Bg.get(c), c.nodeType === 1 && !DA.get(c, "hasDataAttrs"))) {
            for (D = G.length; D--; )
              G[D] && (w = G[D].name, w.indexOf("data-") === 0 && (w = ag(w.slice(5)), yB(c, w, e[w])));
            DA.set(c, "hasDataAttrs", !0);
          }
          return e;
        }
        return typeof Q == "object" ? this.each(function() {
          Bg.set(this, Q);
        }) : HA(this, function(R) {
          var n;
          if (c && R === void 0)
            return n = Bg.get(c, Q), n !== void 0 || (n = yB(c, Q), n !== void 0) ? n : void 0;
          this.each(function() {
            Bg.set(this, Q, R);
          });
        }, null, i, arguments.length > 1, null, !0);
      },
      removeData: function(Q) {
        return this.each(function() {
          Bg.remove(this, Q);
        });
      }
    }), a.extend({
      queue: function(Q, i, D) {
        var w;
        if (Q)
          return i = (i || "fx") + "queue", w = DA.get(Q, i), D && (!w || Array.isArray(D) ? w = DA.access(Q, i, a.makeArray(D)) : w.push(D)), w || [];
      },
      dequeue: function(Q, i) {
        i = i || "fx";
        var D = a.queue(Q, i), w = D.length, e = D.shift(), c = a._queueHooks(Q, i), G = function() {
          a.dequeue(Q, i);
        };
        e === "inprogress" && (e = D.shift(), w--), e && (i === "fx" && D.unshift("inprogress"), delete c.stop, e.call(Q, G, c)), !w && c && c.empty.fire();
      },
      // Not public - generate a queueHooks object, or return the current one
      _queueHooks: function(Q, i) {
        var D = i + "queueHooks";
        return DA.get(Q, D) || DA.access(Q, D, {
          empty: a.Callbacks("once memory").add(function() {
            DA.remove(Q, [i + "queue", D]);
          })
        });
      }
    }), a.fn.extend({
      queue: function(Q, i) {
        var D = 2;
        return typeof Q != "string" && (i = Q, Q = "fx", D--), arguments.length < D ? a.queue(this[0], Q) : i === void 0 ? this : this.each(function() {
          var w = a.queue(this, Q, i);
          a._queueHooks(this, Q), Q === "fx" && w[0] !== "inprogress" && a.dequeue(this, Q);
        });
      },
      dequeue: function(Q) {
        return this.each(function() {
          a.dequeue(this, Q);
        });
      },
      clearQueue: function(Q) {
        return this.queue(Q || "fx", []);
      },
      // Get a promise resolved when queues of a certain type
      // are emptied (fx is the type by default)
      promise: function(Q, i) {
        var D, w = 1, e = a.Deferred(), c = this, G = this.length, R = function() {
          --w || e.resolveWith(c, [c]);
        };
        for (typeof Q != "string" && (i = Q, Q = void 0), Q = Q || "fx"; G--; )
          D = DA.get(c[G], Q + "queueHooks"), D && D.empty && (w++, D.empty.add(R));
        return R(), e.promise(i);
      }
    });
    var rB = /[+-]?(?:\d*\.|)\d+(?:[eE][+-]?\d+|)/.source, QI = new RegExp("^(?:([+-])=|)(" + rB + ")([a-z%]*)$", "i"), Kg = ["Top", "Right", "Bottom", "Left"], dg = q.documentElement, Zg = function(Q) {
      return a.contains(Q.ownerDocument, Q);
    }, dQ = { composed: !0 };
    dg.getRootNode && (Zg = function(Q) {
      return a.contains(Q.ownerDocument, Q) || Q.getRootNode(dQ) === Q.ownerDocument;
    });
    var KI = function(Q, i) {
      return Q = i || Q, Q.style.display === "none" || Q.style.display === "" && // Otherwise, check computed style
      // Support: Firefox <=43 - 45
      // Disconnected elements can have computed display: none, so first confirm that elem is
      // in the document.
      Zg(Q) && a.css(Q, "display") === "none";
    };
    function NB(Q, i, D, w) {
      var e, c, G = 20, R = w ? function() {
        return w.cur();
      } : function() {
        return a.css(Q, i, "");
      }, n = R(), L = D && D[3] || (a.cssNumber[i] ? "" : "px"), p = Q.nodeType && (a.cssNumber[i] || L !== "px" && +n) && QI.exec(a.css(Q, i));
      if (p && p[3] !== L) {
        for (n = n / 2, L = L || p[3], p = +n || 1; G--; )
          a.style(Q, i, p + L), (1 - c) * (1 - (c = R() / n || 0.5)) <= 0 && (G = 0), p = p / c;
        p = p * 2, a.style(Q, i, p + L), D = D || [];
      }
      return D && (p = +p || +n || 0, e = D[1] ? p + (D[1] + 1) * D[2] : +D[2], w && (w.unit = L, w.start = p, w.end = e)), e;
    }
    var kB = {};
    function HQ(Q) {
      var i, D = Q.ownerDocument, w = Q.nodeName, e = kB[w];
      return e || (i = D.body.appendChild(D.createElement(w)), e = a.css(i, "display"), i.parentNode.removeChild(i), e === "none" && (e = "block"), kB[w] = e, e);
    }
    function Tg(Q, i) {
      for (var D, w, e = [], c = 0, G = Q.length; c < G; c++)
        w = Q[c], w.style && (D = w.style.display, i ? (D === "none" && (e[c] = DA.get(w, "display") || null, e[c] || (w.style.display = "")), w.style.display === "" && KI(w) && (e[c] = HQ(w))) : D !== "none" && (e[c] = "none", DA.set(w, "display", D)));
      for (c = 0; c < G; c++)
        e[c] != null && (Q[c].style.display = e[c]);
      return Q;
    }
    a.fn.extend({
      show: function() {
        return Tg(this, !0);
      },
      hide: function() {
        return Tg(this);
      },
      toggle: function(Q) {
        return typeof Q == "boolean" ? Q ? this.show() : this.hide() : this.each(function() {
          KI(this) ? a(this).show() : a(this).hide();
        });
      }
    });
    var CI = /^(?:checkbox|radio)$/i, nB = /<([a-z][^\/\0>\x20\t\r\n\f]*)/i, FB = /^$|^module$|\/(?:java|ecma)script/i;
    (function() {
      var Q = q.createDocumentFragment(), i = Q.appendChild(q.createElement("div")), D = q.createElement("input");
      D.setAttribute("type", "radio"), D.setAttribute("checked", "checked"), D.setAttribute("name", "t"), i.appendChild(D), K.checkClone = i.cloneNode(!0).cloneNode(!0).lastChild.checked, i.innerHTML = "<textarea>x</textarea>", K.noCloneChecked = !!i.cloneNode(!0).lastChild.defaultValue, i.innerHTML = "<option></option>", K.option = !!i.lastChild;
    })();
    var Gg = {
      // XHTML parsers do not magically insert elements in the
      // same way that tag soup parsers do. So we cannot shorten
      // this by omitting <tbody> or other required elements.
      thead: [1, "<table>", "</table>"],
      col: [2, "<table><colgroup>", "</colgroup></table>"],
      tr: [2, "<table><tbody>", "</tbody></table>"],
      td: [3, "<table><tbody><tr>", "</tr></tbody></table>"],
      _default: [0, "", ""]
    };
    Gg.tbody = Gg.tfoot = Gg.colgroup = Gg.caption = Gg.thead, Gg.th = Gg.td, K.option || (Gg.optgroup = Gg.option = [1, "<select multiple='multiple'>", "</select>"]);
    function Qg(Q, i) {
      var D;
      return typeof Q.getElementsByTagName < "u" ? D = Q.getElementsByTagName(i || "*") : typeof Q.querySelectorAll < "u" ? D = Q.querySelectorAll(i || "*") : D = [], i === void 0 || i && yA(Q, i) ? a.merge([Q], D) : D;
    }
    function ZI(Q, i) {
      for (var D = 0, w = Q.length; D < w; D++)
        DA.set(
          Q[D],
          "globalEval",
          !i || DA.get(i[D], "globalEval")
        );
    }
    var pQ = /<|&#?\w+;/;
    function RB(Q, i, D, w, e) {
      for (var c, G, R, n, L, p, V = i.createDocumentFragment(), H = [], u = 0, BA = Q.length; u < BA; u++)
        if (c = Q[u], c || c === 0)
          if (v(c) === "object")
            a.merge(H, c.nodeType ? [c] : c);
          else if (!pQ.test(c))
            H.push(i.createTextNode(c));
          else {
            for (G = G || V.appendChild(i.createElement("div")), R = (nB.exec(c) || ["", ""])[1].toLowerCase(), n = Gg[R] || Gg._default, G.innerHTML = n[1] + a.htmlPrefilter(c) + n[2], p = n[0]; p--; )
              G = G.lastChild;
            a.merge(H, G.childNodes), G = V.firstChild, G.textContent = "";
          }
      for (V.textContent = "", u = 0; c = H[u++]; ) {
        if (w && a.inArray(c, w) > -1) {
          e && e.push(c);
          continue;
        }
        if (L = Zg(c), G = Qg(V.appendChild(c), "script"), L && ZI(G), D)
          for (p = 0; c = G[p++]; )
            FB.test(c.type || "") && D.push(c);
      }
      return V;
    }
    var KB = /^([^.]*)(?:\.(.+)|)/;
    function Wg() {
      return !0;
    }
    function Vg() {
      return !1;
    }
    function fQ(Q, i) {
      return Q === uQ() == (i === "focus");
    }
    function uQ() {
      try {
        return q.activeElement;
      } catch {
      }
    }
    function TI(Q, i, D, w, e, c) {
      var G, R;
      if (typeof i == "object") {
        typeof D != "string" && (w = w || D, D = void 0);
        for (R in i)
          TI(Q, R, D, w, i[R], c);
        return Q;
      }
      if (w == null && e == null ? (e = D, w = D = void 0) : e == null && (typeof D == "string" ? (e = w, w = void 0) : (e = w, w = D, D = void 0)), e === !1)
        e = Vg;
      else if (!e)
        return Q;
      return c === 1 && (G = e, e = function(n) {
        return a().off(n), G.apply(this, arguments);
      }, e.guid = G.guid || (G.guid = a.guid++)), Q.each(function() {
        a.event.add(this, i, e, w, D);
      });
    }
    a.event = {
      global: {},
      add: function(Q, i, D, w, e) {
        var c, G, R, n, L, p, V, H, u, BA, cA, QA = DA.get(Q);
        if (II(Q))
          for (D.handler && (c = D, D = c.handler, e = c.selector), e && a.find.matchesSelector(dg, e), D.guid || (D.guid = a.guid++), (n = QA.events) || (n = QA.events = /* @__PURE__ */ Object.create(null)), (G = QA.handle) || (G = QA.handle = function(TA) {
            return typeof a < "u" && a.event.triggered !== TA.type ? a.event.dispatch.apply(Q, arguments) : void 0;
          }), i = (i || "").match(Dg) || [""], L = i.length; L--; )
            R = KB.exec(i[L]) || [], u = cA = R[1], BA = (R[2] || "").split(".").sort(), u && (V = a.event.special[u] || {}, u = (e ? V.delegateType : V.bindType) || u, V = a.event.special[u] || {}, p = a.extend({
              type: u,
              origType: cA,
              data: w,
              handler: D,
              guid: D.guid,
              selector: e,
              needsContext: e && a.expr.match.needsContext.test(e),
              namespace: BA.join(".")
            }, c), (H = n[u]) || (H = n[u] = [], H.delegateCount = 0, (!V.setup || V.setup.call(Q, w, BA, G) === !1) && Q.addEventListener && Q.addEventListener(u, G)), V.add && (V.add.call(Q, p), p.handler.guid || (p.handler.guid = D.guid)), e ? H.splice(H.delegateCount++, 0, p) : H.push(p), a.event.global[u] = !0);
      },
      // Detach an event or set of events from an element
      remove: function(Q, i, D, w, e) {
        var c, G, R, n, L, p, V, H, u, BA, cA, QA = DA.hasData(Q) && DA.get(Q);
        if (!(!QA || !(n = QA.events))) {
          for (i = (i || "").match(Dg) || [""], L = i.length; L--; ) {
            if (R = KB.exec(i[L]) || [], u = cA = R[1], BA = (R[2] || "").split(".").sort(), !u) {
              for (u in n)
                a.event.remove(Q, u + i[L], D, w, !0);
              continue;
            }
            for (V = a.event.special[u] || {}, u = (w ? V.delegateType : V.bindType) || u, H = n[u] || [], R = R[2] && new RegExp("(^|\\.)" + BA.join("\\.(?:.*\\.|)") + "(\\.|$)"), G = c = H.length; c--; )
              p = H[c], (e || cA === p.origType) && (!D || D.guid === p.guid) && (!R || R.test(p.namespace)) && (!w || w === p.selector || w === "**" && p.selector) && (H.splice(c, 1), p.selector && H.delegateCount--, V.remove && V.remove.call(Q, p));
            G && !H.length && ((!V.teardown || V.teardown.call(Q, BA, QA.handle) === !1) && a.removeEvent(Q, u, QA.handle), delete n[u]);
          }
          a.isEmptyObject(n) && DA.remove(Q, "handle events");
        }
      },
      dispatch: function(Q) {
        var i, D, w, e, c, G, R = new Array(arguments.length), n = a.event.fix(Q), L = (DA.get(this, "events") || /* @__PURE__ */ Object.create(null))[n.type] || [], p = a.event.special[n.type] || {};
        for (R[0] = n, i = 1; i < arguments.length; i++)
          R[i] = arguments[i];
        if (n.delegateTarget = this, !(p.preDispatch && p.preDispatch.call(this, n) === !1)) {
          for (G = a.event.handlers.call(this, n, L), i = 0; (e = G[i++]) && !n.isPropagationStopped(); )
            for (n.currentTarget = e.elem, D = 0; (c = e.handlers[D++]) && !n.isImmediatePropagationStopped(); )
              (!n.rnamespace || c.namespace === !1 || n.rnamespace.test(c.namespace)) && (n.handleObj = c, n.data = c.data, w = ((a.event.special[c.origType] || {}).handle || c.handler).apply(e.elem, R), w !== void 0 && (n.result = w) === !1 && (n.preventDefault(), n.stopPropagation()));
          return p.postDispatch && p.postDispatch.call(this, n), n.result;
        }
      },
      handlers: function(Q, i) {
        var D, w, e, c, G, R = [], n = i.delegateCount, L = Q.target;
        if (n && // Support: IE <=9
        // Black-hole SVG <use> instance trees (trac-13180)
        L.nodeType && // Support: Firefox <=42
        // Suppress spec-violating clicks indicating a non-primary pointer button (trac-3861)
        // https://www.w3.org/TR/DOM-Level-3-Events/#event-type-click
        // Support: IE 11 only
        // ...but not arrow key "clicks" of radio inputs, which can have `button` -1 (gh-2343)
        !(Q.type === "click" && Q.button >= 1)) {
          for (; L !== this; L = L.parentNode || this)
            if (L.nodeType === 1 && !(Q.type === "click" && L.disabled === !0)) {
              for (c = [], G = {}, D = 0; D < n; D++)
                w = i[D], e = w.selector + " ", G[e] === void 0 && (G[e] = w.needsContext ? a(e, this).index(L) > -1 : a.find(e, this, null, [L]).length), G[e] && c.push(w);
              c.length && R.push({ elem: L, handlers: c });
            }
        }
        return L = this, n < i.length && R.push({ elem: L, handlers: i.slice(n) }), R;
      },
      addProp: function(Q, i) {
        Object.defineProperty(a.Event.prototype, Q, {
          enumerable: !0,
          configurable: !0,
          get: l(i) ? function() {
            if (this.originalEvent)
              return i(this.originalEvent);
          } : function() {
            if (this.originalEvent)
              return this.originalEvent[Q];
          },
          set: function(D) {
            Object.defineProperty(this, Q, {
              enumerable: !0,
              configurable: !0,
              writable: !0,
              value: D
            });
          }
        });
      },
      fix: function(Q) {
        return Q[a.expando] ? Q : new a.Event(Q);
      },
      special: {
        load: {
          // Prevent triggered image.load events from bubbling to window.load
          noBubble: !0
        },
        click: {
          // Utilize native event to ensure correct state for checkable inputs
          setup: function(Q) {
            var i = this || Q;
            return CI.test(i.type) && i.click && yA(i, "input") && YI(i, "click", Wg), !1;
          },
          trigger: function(Q) {
            var i = this || Q;
            return CI.test(i.type) && i.click && yA(i, "input") && YI(i, "click"), !0;
          },
          // For cross-browser consistency, suppress native .click() on links
          // Also prevent it if we're currently inside a leveraged native-event stack
          _default: function(Q) {
            var i = Q.target;
            return CI.test(i.type) && i.click && yA(i, "input") && DA.get(i, "click") || yA(i, "a");
          }
        },
        beforeunload: {
          postDispatch: function(Q) {
            Q.result !== void 0 && Q.originalEvent && (Q.originalEvent.returnValue = Q.result);
          }
        }
      }
    };
    function YI(Q, i, D) {
      if (!D) {
        DA.get(Q, i) === void 0 && a.event.add(Q, i, Wg);
        return;
      }
      DA.set(Q, i, !1), a.event.add(Q, i, {
        namespace: !1,
        handler: function(w) {
          var e, c, G = DA.get(this, i);
          if (w.isTrigger & 1 && this[i]) {
            if (G.length)
              (a.event.special[i] || {}).delegateType && w.stopPropagation();
            else if (G = E.call(arguments), DA.set(this, i, G), e = D(this, i), this[i](), c = DA.get(this, i), G !== c || e ? DA.set(this, i, !1) : c = {}, G !== c)
              return w.stopImmediatePropagation(), w.preventDefault(), c && c.value;
          } else
            G.length && (DA.set(this, i, {
              value: a.event.trigger(
                // Support: IE <=9 - 11+
                // Extend with the prototype to reset the above stopImmediatePropagation()
                a.extend(G[0], a.Event.prototype),
                G.slice(1),
                this
              )
            }), w.stopImmediatePropagation());
        }
      });
    }
    a.removeEvent = function(Q, i, D) {
      Q.removeEventListener && Q.removeEventListener(i, D);
    }, a.Event = function(Q, i) {
      if (!(this instanceof a.Event))
        return new a.Event(Q, i);
      Q && Q.type ? (this.originalEvent = Q, this.type = Q.type, this.isDefaultPrevented = Q.defaultPrevented || Q.defaultPrevented === void 0 && // Support: Android <=2.3 only
      Q.returnValue === !1 ? Wg : Vg, this.target = Q.target && Q.target.nodeType === 3 ? Q.target.parentNode : Q.target, this.currentTarget = Q.currentTarget, this.relatedTarget = Q.relatedTarget) : this.type = Q, i && a.extend(this, i), this.timeStamp = Q && Q.timeStamp || Date.now(), this[a.expando] = !0;
    }, a.Event.prototype = {
      constructor: a.Event,
      isDefaultPrevented: Vg,
      isPropagationStopped: Vg,
      isImmediatePropagationStopped: Vg,
      isSimulated: !1,
      preventDefault: function() {
        var Q = this.originalEvent;
        this.isDefaultPrevented = Wg, Q && !this.isSimulated && Q.preventDefault();
      },
      stopPropagation: function() {
        var Q = this.originalEvent;
        this.isPropagationStopped = Wg, Q && !this.isSimulated && Q.stopPropagation();
      },
      stopImmediatePropagation: function() {
        var Q = this.originalEvent;
        this.isImmediatePropagationStopped = Wg, Q && !this.isSimulated && Q.stopImmediatePropagation(), this.stopPropagation();
      }
    }, a.each({
      altKey: !0,
      bubbles: !0,
      cancelable: !0,
      changedTouches: !0,
      ctrlKey: !0,
      detail: !0,
      eventPhase: !0,
      metaKey: !0,
      pageX: !0,
      pageY: !0,
      shiftKey: !0,
      view: !0,
      char: !0,
      code: !0,
      charCode: !0,
      key: !0,
      keyCode: !0,
      button: !0,
      buttons: !0,
      clientX: !0,
      clientY: !0,
      offsetX: !0,
      offsetY: !0,
      pointerId: !0,
      pointerType: !0,
      screenX: !0,
      screenY: !0,
      targetTouches: !0,
      toElement: !0,
      touches: !0,
      which: !0
    }, a.event.addProp), a.each({ focus: "focusin", blur: "focusout" }, function(Q, i) {
      a.event.special[Q] = {
        // Utilize native event if possible so blur/focus sequence is correct
        setup: function() {
          return YI(this, Q, fQ), !1;
        },
        trigger: function() {
          return YI(this, Q), !0;
        },
        // Suppress native focus or blur if we're currently inside
        // a leveraged native-event stack
        _default: function(D) {
          return DA.get(D.target, Q);
        },
        delegateType: i
      };
    }), a.each({
      mouseenter: "mouseover",
      mouseleave: "mouseout",
      pointerenter: "pointerover",
      pointerleave: "pointerout"
    }, function(Q, i) {
      a.event.special[Q] = {
        delegateType: i,
        bindType: i,
        handle: function(D) {
          var w, e = this, c = D.relatedTarget, G = D.handleObj;
          return (!c || c !== e && !a.contains(e, c)) && (D.type = G.origType, w = G.handler.apply(this, arguments), D.type = i), w;
        }
      };
    }), a.fn.extend({
      on: function(Q, i, D, w) {
        return TI(this, Q, i, D, w);
      },
      one: function(Q, i, D, w) {
        return TI(this, Q, i, D, w, 1);
      },
      off: function(Q, i, D) {
        var w, e;
        if (Q && Q.preventDefault && Q.handleObj)
          return w = Q.handleObj, a(Q.delegateTarget).off(
            w.namespace ? w.origType + "." + w.namespace : w.origType,
            w.selector,
            w.handler
          ), this;
        if (typeof Q == "object") {
          for (e in Q)
            this.off(e, i, Q[e]);
          return this;
        }
        return (i === !1 || typeof i == "function") && (D = i, i = void 0), D === !1 && (D = Vg), this.each(function() {
          a.event.remove(this, Q, D, i);
        });
      }
    });
    var xQ = /<script|<style|<link/i, mQ = /checked\s*(?:[^=]|=\s*.checked.)/i, OQ = /^\s*<!\[CDATA\[|\]\]>\s*$/g;
    function YB(Q, i) {
      return yA(Q, "table") && yA(i.nodeType !== 11 ? i : i.firstChild, "tr") && a(Q).children("tbody")[0] || Q;
    }
    function vQ(Q) {
      return Q.type = (Q.getAttribute("type") !== null) + "/" + Q.type, Q;
    }
    function bQ(Q) {
      return (Q.type || "").slice(0, 5) === "true/" ? Q.type = Q.type.slice(5) : Q.removeAttribute("type"), Q;
    }
    function lB(Q, i) {
      var D, w, e, c, G, R, n;
      if (i.nodeType === 1) {
        if (DA.hasData(Q) && (c = DA.get(Q), n = c.events, n)) {
          DA.remove(i, "handle events");
          for (e in n)
            for (D = 0, w = n[e].length; D < w; D++)
              a.event.add(i, e, n[e][D]);
        }
        Bg.hasData(Q) && (G = Bg.access(Q), R = a.extend({}, G), Bg.set(i, R));
      }
    }
    function ZQ(Q, i) {
      var D = i.nodeName.toLowerCase();
      D === "input" && CI.test(Q.type) ? i.checked = Q.checked : (D === "input" || D === "textarea") && (i.defaultValue = Q.defaultValue);
    }
    function jg(Q, i, D, w) {
      i = o(i);
      var e, c, G, R, n, L, p = 0, V = Q.length, H = V - 1, u = i[0], BA = l(u);
      if (BA || V > 1 && typeof u == "string" && !K.checkClone && mQ.test(u))
        return Q.each(function(cA) {
          var QA = Q.eq(cA);
          BA && (i[0] = u.call(this, cA, QA.html())), jg(QA, i, D, w);
        });
      if (V && (e = RB(i, Q[0].ownerDocument, !1, Q, w), c = e.firstChild, e.childNodes.length === 1 && (e = c), c || w)) {
        for (G = a.map(Qg(e, "script"), vQ), R = G.length; p < V; p++)
          n = e, p !== H && (n = a.clone(n, !0, !0), R && a.merge(G, Qg(n, "script"))), D.call(Q[p], n, p);
        if (R)
          for (L = G[G.length - 1].ownerDocument, a.map(G, bQ), p = 0; p < R; p++)
            n = G[p], FB.test(n.type || "") && !DA.access(n, "globalEval") && a.contains(L, n) && (n.src && (n.type || "").toLowerCase() !== "module" ? a._evalUrl && !n.noModule && a._evalUrl(n.src, {
              nonce: n.nonce || n.getAttribute("nonce")
            }, L) : T(n.textContent.replace(OQ, ""), n, L));
      }
      return Q;
    }
    function JB(Q, i, D) {
      for (var w, e = i ? a.filter(i, Q) : Q, c = 0; (w = e[c]) != null; c++)
        !D && w.nodeType === 1 && a.cleanData(Qg(w)), w.parentNode && (D && Zg(w) && ZI(Qg(w, "script")), w.parentNode.removeChild(w));
      return Q;
    }
    a.extend({
      htmlPrefilter: function(Q) {
        return Q;
      },
      clone: function(Q, i, D) {
        var w, e, c, G, R = Q.cloneNode(!0), n = Zg(Q);
        if (!K.noCloneChecked && (Q.nodeType === 1 || Q.nodeType === 11) && !a.isXMLDoc(Q))
          for (G = Qg(R), c = Qg(Q), w = 0, e = c.length; w < e; w++)
            ZQ(c[w], G[w]);
        if (i)
          if (D)
            for (c = c || Qg(Q), G = G || Qg(R), w = 0, e = c.length; w < e; w++)
              lB(c[w], G[w]);
          else
            lB(Q, R);
        return G = Qg(R, "script"), G.length > 0 && ZI(G, !n && Qg(Q, "script")), R;
      },
      cleanData: function(Q) {
        for (var i, D, w, e = a.event.special, c = 0; (D = Q[c]) !== void 0; c++)
          if (II(D)) {
            if (i = D[DA.expando]) {
              if (i.events)
                for (w in i.events)
                  e[w] ? a.event.remove(D, w) : a.removeEvent(D, w, i.handle);
              D[DA.expando] = void 0;
            }
            D[Bg.expando] && (D[Bg.expando] = void 0);
          }
      }
    }), a.fn.extend({
      detach: function(Q) {
        return JB(this, Q, !0);
      },
      remove: function(Q) {
        return JB(this, Q);
      },
      text: function(Q) {
        return HA(this, function(i) {
          return i === void 0 ? a.text(this) : this.empty().each(function() {
            (this.nodeType === 1 || this.nodeType === 11 || this.nodeType === 9) && (this.textContent = i);
          });
        }, null, Q, arguments.length);
      },
      append: function() {
        return jg(this, arguments, function(Q) {
          if (this.nodeType === 1 || this.nodeType === 11 || this.nodeType === 9) {
            var i = YB(this, Q);
            i.appendChild(Q);
          }
        });
      },
      prepend: function() {
        return jg(this, arguments, function(Q) {
          if (this.nodeType === 1 || this.nodeType === 11 || this.nodeType === 9) {
            var i = YB(this, Q);
            i.insertBefore(Q, i.firstChild);
          }
        });
      },
      before: function() {
        return jg(this, arguments, function(Q) {
          this.parentNode && this.parentNode.insertBefore(Q, this);
        });
      },
      after: function() {
        return jg(this, arguments, function(Q) {
          this.parentNode && this.parentNode.insertBefore(Q, this.nextSibling);
        });
      },
      empty: function() {
        for (var Q, i = 0; (Q = this[i]) != null; i++)
          Q.nodeType === 1 && (a.cleanData(Qg(Q, !1)), Q.textContent = "");
        return this;
      },
      clone: function(Q, i) {
        return Q = Q ?? !1, i = i ?? Q, this.map(function() {
          return a.clone(this, Q, i);
        });
      },
      html: function(Q) {
        return HA(this, function(i) {
          var D = this[0] || {}, w = 0, e = this.length;
          if (i === void 0 && D.nodeType === 1)
            return D.innerHTML;
          if (typeof i == "string" && !xQ.test(i) && !Gg[(nB.exec(i) || ["", ""])[1].toLowerCase()]) {
            i = a.htmlPrefilter(i);
            try {
              for (; w < e; w++)
                D = this[w] || {}, D.nodeType === 1 && (a.cleanData(Qg(D, !1)), D.innerHTML = i);
              D = 0;
            } catch {
            }
          }
          D && this.empty().append(i);
        }, null, Q, arguments.length);
      },
      replaceWith: function() {
        var Q = [];
        return jg(this, arguments, function(i) {
          var D = this.parentNode;
          a.inArray(this, Q) < 0 && (a.cleanData(Qg(this)), D && D.replaceChild(i, this));
        }, Q);
      }
    }), a.each({
      appendTo: "append",
      prependTo: "prepend",
      insertBefore: "before",
      insertAfter: "after",
      replaceAll: "replaceWith"
    }, function(Q, i) {
      a.fn[Q] = function(D) {
        for (var w, e = [], c = a(D), G = c.length - 1, R = 0; R <= G; R++)
          w = R === G ? this : this.clone(!0), a(c[R])[i](w), t.apply(e, w.get());
        return this.pushStack(e);
      };
    });
    var WI = new RegExp("^(" + rB + ")(?!px)[a-z%]+$", "i"), VI = /^--/, lI = function(Q) {
      var i = Q.ownerDocument.defaultView;
      return (!i || !i.opener) && (i = A), i.getComputedStyle(Q);
    }, UB = function(Q, i, D) {
      var w, e, c = {};
      for (e in i)
        c[e] = Q.style[e], Q.style[e] = i[e];
      w = D.call(Q);
      for (e in i)
        Q.style[e] = c[e];
      return w;
    }, TQ = new RegExp(Kg.join("|"), "i"), SB = "[\\x20\\t\\r\\n\\f]", WQ = new RegExp(
      "^" + SB + "+|((?:^|[^\\\\])(?:\\\\.)*)" + SB + "+$",
      "g"
    );
    (function() {
      function Q() {
        if (L) {
          n.style.cssText = "position:absolute;left:-11111px;width:60px;margin-top:1px;padding:0;border:0", L.style.cssText = "position:relative;display:block;box-sizing:border-box;overflow:scroll;margin:auto;border:1px;padding:1px;width:60%;top:1%", dg.appendChild(n).appendChild(L);
          var p = A.getComputedStyle(L);
          D = p.top !== "1%", R = i(p.marginLeft) === 12, L.style.right = "60%", c = i(p.right) === 36, w = i(p.width) === 36, L.style.position = "absolute", e = i(L.offsetWidth / 3) === 12, dg.removeChild(n), L = null;
        }
      }
      function i(p) {
        return Math.round(parseFloat(p));
      }
      var D, w, e, c, G, R, n = q.createElement("div"), L = q.createElement("div");
      L.style && (L.style.backgroundClip = "content-box", L.cloneNode(!0).style.backgroundClip = "", K.clearCloneStyle = L.style.backgroundClip === "content-box", a.extend(K, {
        boxSizingReliable: function() {
          return Q(), w;
        },
        pixelBoxStyles: function() {
          return Q(), c;
        },
        pixelPosition: function() {
          return Q(), D;
        },
        reliableMarginLeft: function() {
          return Q(), R;
        },
        scrollboxSize: function() {
          return Q(), e;
        },
        // Support: IE 9 - 11+, Edge 15 - 18+
        // IE/Edge misreport `getComputedStyle` of table rows with width/height
        // set in CSS while `offset*` properties report correct values.
        // Behavior in IE 9 is more subtle than in newer versions & it passes
        // some versions of this test; make sure not to make it pass there!
        //
        // Support: Firefox 70+
        // Only Firefox includes border widths
        // in computed dimensions. (gh-4529)
        reliableTrDimensions: function() {
          var p, V, H, u;
          return G == null && (p = q.createElement("table"), V = q.createElement("tr"), H = q.createElement("div"), p.style.cssText = "position:absolute;left:-11111px;border-collapse:separate", V.style.cssText = "border:1px solid", V.style.height = "1px", H.style.height = "9px", H.style.display = "block", dg.appendChild(p).appendChild(V).appendChild(H), u = A.getComputedStyle(V), G = parseInt(u.height, 10) + parseInt(u.borderTopWidth, 10) + parseInt(u.borderBottomWidth, 10) === V.offsetHeight, dg.removeChild(p)), G;
        }
      }));
    })();
    function EI(Q, i, D) {
      var w, e, c, G, R = VI.test(i), n = Q.style;
      return D = D || lI(Q), D && (G = D.getPropertyValue(i) || D[i], R && (G = G.replace(WQ, "$1")), G === "" && !Zg(Q) && (G = a.style(Q, i)), !K.pixelBoxStyles() && WI.test(G) && TQ.test(i) && (w = n.width, e = n.minWidth, c = n.maxWidth, n.minWidth = n.maxWidth = n.width = G, G = D.width, n.width = w, n.minWidth = e, n.maxWidth = c)), G !== void 0 ? (
        // Support: IE <=9 - 11 only
        // IE returns zIndex value as an integer.
        G + ""
      ) : G;
    }
    function LB(Q, i) {
      return {
        get: function() {
          if (Q()) {
            delete this.get;
            return;
          }
          return (this.get = i).apply(this, arguments);
        }
      };
    }
    var qB = ["Webkit", "Moz", "ms"], dB = q.createElement("div").style, HB = {};
    function VQ(Q) {
      for (var i = Q[0].toUpperCase() + Q.slice(1), D = qB.length; D--; )
        if (Q = qB[D] + i, Q in dB)
          return Q;
    }
    function jI(Q) {
      var i = a.cssProps[Q] || HB[Q];
      return i || (Q in dB ? Q : HB[Q] = VQ(Q) || Q);
    }
    var jQ = /^(none|table(?!-c[ea]).+)/, _Q = { position: "absolute", visibility: "hidden", display: "block" }, pB = {
      letterSpacing: "0",
      fontWeight: "400"
    };
    function fB(Q, i, D) {
      var w = QI.exec(i);
      return w ? (
        // Guard against undefined "subtract", e.g., when used as in cssHooks
        Math.max(0, w[2] - (D || 0)) + (w[3] || "px")
      ) : i;
    }
    function _I(Q, i, D, w, e, c) {
      var G = i === "width" ? 1 : 0, R = 0, n = 0;
      if (D === (w ? "border" : "content"))
        return 0;
      for (; G < 4; G += 2)
        D === "margin" && (n += a.css(Q, D + Kg[G], !0, e)), w ? (D === "content" && (n -= a.css(Q, "padding" + Kg[G], !0, e)), D !== "margin" && (n -= a.css(Q, "border" + Kg[G] + "Width", !0, e))) : (n += a.css(Q, "padding" + Kg[G], !0, e), D !== "padding" ? n += a.css(Q, "border" + Kg[G] + "Width", !0, e) : R += a.css(Q, "border" + Kg[G] + "Width", !0, e));
      return !w && c >= 0 && (n += Math.max(0, Math.ceil(
        Q["offset" + i[0].toUpperCase() + i.slice(1)] - c - n - R - 0.5
        // If offsetWidth/offsetHeight is unknown, then we can't determine content-box scroll gutter
        // Use an explicit zero to avoid NaN (gh-3964)
      )) || 0), n;
    }
    function uB(Q, i, D) {
      var w = lI(Q), e = !K.boxSizingReliable() || D, c = e && a.css(Q, "boxSizing", !1, w) === "border-box", G = c, R = EI(Q, i, w), n = "offset" + i[0].toUpperCase() + i.slice(1);
      if (WI.test(R)) {
        if (!D)
          return R;
        R = "auto";
      }
      return (!K.boxSizingReliable() && c || // Support: IE 10 - 11+, Edge 15 - 18+
      // IE/Edge misreport `getComputedStyle` of table rows with width/height
      // set in CSS while `offset*` properties report correct values.
      // Interestingly, in some cases IE 9 doesn't suffer from this issue.
      !K.reliableTrDimensions() && yA(Q, "tr") || // Fall back to offsetWidth/offsetHeight when value is "auto"
      // This happens for inline elements with no explicit setting (gh-3571)
      R === "auto" || // Support: Android <=4.1 - 4.3 only
      // Also use offsetWidth/offsetHeight for misreported inline dimensions (gh-3602)
      !parseFloat(R) && a.css(Q, "display", !1, w) === "inline") && // Make sure the element is visible & connected
      Q.getClientRects().length && (c = a.css(Q, "boxSizing", !1, w) === "border-box", G = n in Q, G && (R = Q[n])), R = parseFloat(R) || 0, R + _I(
        Q,
        i,
        D || (c ? "border" : "content"),
        G,
        w,
        // Provide the current computed size to request scroll gutter calculation (gh-3589)
        R
      ) + "px";
    }
    a.extend({
      // Add in style property hooks for overriding the default
      // behavior of getting and setting a style property
      cssHooks: {
        opacity: {
          get: function(Q, i) {
            if (i) {
              var D = EI(Q, "opacity");
              return D === "" ? "1" : D;
            }
          }
        }
      },
      // Don't automatically add "px" to these possibly-unitless properties
      cssNumber: {
        animationIterationCount: !0,
        columnCount: !0,
        fillOpacity: !0,
        flexGrow: !0,
        flexShrink: !0,
        fontWeight: !0,
        gridArea: !0,
        gridColumn: !0,
        gridColumnEnd: !0,
        gridColumnStart: !0,
        gridRow: !0,
        gridRowEnd: !0,
        gridRowStart: !0,
        lineHeight: !0,
        opacity: !0,
        order: !0,
        orphans: !0,
        widows: !0,
        zIndex: !0,
        zoom: !0
      },
      // Add in properties whose names you wish to fix before
      // setting or getting the value
      cssProps: {},
      // Get and set the style property on a DOM Node
      style: function(Q, i, D, w) {
        if (!(!Q || Q.nodeType === 3 || Q.nodeType === 8 || !Q.style)) {
          var e, c, G, R = ag(i), n = VI.test(i), L = Q.style;
          if (n || (i = jI(R)), G = a.cssHooks[i] || a.cssHooks[R], D !== void 0) {
            if (c = typeof D, c === "string" && (e = QI.exec(D)) && e[1] && (D = NB(Q, i, e), c = "number"), D == null || D !== D)
              return;
            c === "number" && !n && (D += e && e[3] || (a.cssNumber[R] ? "" : "px")), !K.clearCloneStyle && D === "" && i.indexOf("background") === 0 && (L[i] = "inherit"), (!G || !("set" in G) || (D = G.set(Q, D, w)) !== void 0) && (n ? L.setProperty(i, D) : L[i] = D);
          } else
            return G && "get" in G && (e = G.get(Q, !1, w)) !== void 0 ? e : L[i];
        }
      },
      css: function(Q, i, D, w) {
        var e, c, G, R = ag(i), n = VI.test(i);
        return n || (i = jI(R)), G = a.cssHooks[i] || a.cssHooks[R], G && "get" in G && (e = G.get(Q, !0, D)), e === void 0 && (e = EI(Q, i, w)), e === "normal" && i in pB && (e = pB[i]), D === "" || D ? (c = parseFloat(e), D === !0 || isFinite(c) ? c || 0 : e) : e;
      }
    }), a.each(["height", "width"], function(Q, i) {
      a.cssHooks[i] = {
        get: function(D, w, e) {
          if (w)
            return jQ.test(a.css(D, "display")) && // Support: Safari 8+
            // Table columns in Safari have non-zero offsetWidth & zero
            // getBoundingClientRect().width unless display is changed.
            // Support: IE <=11 only
            // Running getBoundingClientRect on a disconnected node
            // in IE throws an error.
            (!D.getClientRects().length || !D.getBoundingClientRect().width) ? UB(D, _Q, function() {
              return uB(D, i, e);
            }) : uB(D, i, e);
        },
        set: function(D, w, e) {
          var c, G = lI(D), R = !K.scrollboxSize() && G.position === "absolute", n = R || e, L = n && a.css(D, "boxSizing", !1, G) === "border-box", p = e ? _I(
            D,
            i,
            e,
            L,
            G
          ) : 0;
          return L && R && (p -= Math.ceil(
            D["offset" + i[0].toUpperCase() + i.slice(1)] - parseFloat(G[i]) - _I(D, i, "border", !1, G) - 0.5
          )), p && (c = QI.exec(w)) && (c[3] || "px") !== "px" && (D.style[i] = w, w = a.css(D, i)), fB(D, w, p);
        }
      };
    }), a.cssHooks.marginLeft = LB(
      K.reliableMarginLeft,
      function(Q, i) {
        if (i)
          return (parseFloat(EI(Q, "marginLeft")) || Q.getBoundingClientRect().left - UB(Q, { marginLeft: 0 }, function() {
            return Q.getBoundingClientRect().left;
          })) + "px";
      }
    ), a.each({
      margin: "",
      padding: "",
      border: "Width"
    }, function(Q, i) {
      a.cssHooks[Q + i] = {
        expand: function(D) {
          for (var w = 0, e = {}, c = typeof D == "string" ? D.split(" ") : [D]; w < 4; w++)
            e[Q + Kg[w] + i] = c[w] || c[w - 2] || c[0];
          return e;
        }
      }, Q !== "margin" && (a.cssHooks[Q + i].set = fB);
    }), a.fn.extend({
      css: function(Q, i) {
        return HA(this, function(D, w, e) {
          var c, G, R = {}, n = 0;
          if (Array.isArray(w)) {
            for (c = lI(D), G = w.length; n < G; n++)
              R[w[n]] = a.css(D, w[n], !1, c);
            return R;
          }
          return e !== void 0 ? a.style(D, w, e) : a.css(D, w);
        }, Q, i, arguments.length > 1);
      }
    });
    function Cg(Q, i, D, w, e) {
      return new Cg.prototype.init(Q, i, D, w, e);
    }
    a.Tween = Cg, Cg.prototype = {
      constructor: Cg,
      init: function(Q, i, D, w, e, c) {
        this.elem = Q, this.prop = D, this.easing = e || a.easing._default, this.options = i, this.start = this.now = this.cur(), this.end = w, this.unit = c || (a.cssNumber[D] ? "" : "px");
      },
      cur: function() {
        var Q = Cg.propHooks[this.prop];
        return Q && Q.get ? Q.get(this) : Cg.propHooks._default.get(this);
      },
      run: function(Q) {
        var i, D = Cg.propHooks[this.prop];
        return this.options.duration ? this.pos = i = a.easing[this.easing](
          Q,
          this.options.duration * Q,
          0,
          1,
          this.options.duration
        ) : this.pos = i = Q, this.now = (this.end - this.start) * i + this.start, this.options.step && this.options.step.call(this.elem, this.now, this), D && D.set ? D.set(this) : Cg.propHooks._default.set(this), this;
      }
    }, Cg.prototype.init.prototype = Cg.prototype, Cg.propHooks = {
      _default: {
        get: function(Q) {
          var i;
          return Q.elem.nodeType !== 1 || Q.elem[Q.prop] != null && Q.elem.style[Q.prop] == null ? Q.elem[Q.prop] : (i = a.css(Q.elem, Q.prop, ""), !i || i === "auto" ? 0 : i);
        },
        set: function(Q) {
          a.fx.step[Q.prop] ? a.fx.step[Q.prop](Q) : Q.elem.nodeType === 1 && (a.cssHooks[Q.prop] || Q.elem.style[jI(Q.prop)] != null) ? a.style(Q.elem, Q.prop, Q.now + Q.unit) : Q.elem[Q.prop] = Q.now;
        }
      }
    }, Cg.propHooks.scrollTop = Cg.propHooks.scrollLeft = {
      set: function(Q) {
        Q.elem.nodeType && Q.elem.parentNode && (Q.elem[Q.prop] = Q.now);
      }
    }, a.easing = {
      linear: function(Q) {
        return Q;
      },
      swing: function(Q) {
        return 0.5 - Math.cos(Q * Math.PI) / 2;
      },
      _default: "swing"
    }, a.fx = Cg.prototype.init, a.fx.step = {};
    var _g, JI, PQ = /^(?:toggle|show|hide)$/, XQ = /queueHooks$/;
    function PI() {
      JI && (q.hidden === !1 && A.requestAnimationFrame ? A.requestAnimationFrame(PI) : A.setTimeout(PI, a.fx.interval), a.fx.tick());
    }
    function xB() {
      return A.setTimeout(function() {
        _g = void 0;
      }), _g = Date.now();
    }
    function UI(Q, i) {
      var D, w = 0, e = { height: Q };
      for (i = i ? 1 : 0; w < 4; w += 2 - i)
        D = Kg[w], e["margin" + D] = e["padding" + D] = Q;
      return i && (e.opacity = e.width = Q), e;
    }
    function mB(Q, i, D) {
      for (var w, e = (yg.tweeners[i] || []).concat(yg.tweeners["*"]), c = 0, G = e.length; c < G; c++)
        if (w = e[c].call(D, i, Q))
          return w;
    }
    function zQ(Q, i, D) {
      var w, e, c, G, R, n, L, p, V = "width" in i || "height" in i, H = this, u = {}, BA = Q.style, cA = Q.nodeType && KI(Q), QA = DA.get(Q, "fxshow");
      D.queue || (G = a._queueHooks(Q, "fx"), G.unqueued == null && (G.unqueued = 0, R = G.empty.fire, G.empty.fire = function() {
        G.unqueued || R();
      }), G.unqueued++, H.always(function() {
        H.always(function() {
          G.unqueued--, a.queue(Q, "fx").length || G.empty.fire();
        });
      }));
      for (w in i)
        if (e = i[w], PQ.test(e)) {
          if (delete i[w], c = c || e === "toggle", e === (cA ? "hide" : "show"))
            if (e === "show" && QA && QA[w] !== void 0)
              cA = !0;
            else
              continue;
          u[w] = QA && QA[w] || a.style(Q, w);
        }
      if (n = !a.isEmptyObject(i), !(!n && a.isEmptyObject(u))) {
        V && Q.nodeType === 1 && (D.overflow = [BA.overflow, BA.overflowX, BA.overflowY], L = QA && QA.display, L == null && (L = DA.get(Q, "display")), p = a.css(Q, "display"), p === "none" && (L ? p = L : (Tg([Q], !0), L = Q.style.display || L, p = a.css(Q, "display"), Tg([Q]))), (p === "inline" || p === "inline-block" && L != null) && a.css(Q, "float") === "none" && (n || (H.done(function() {
          BA.display = L;
        }), L == null && (p = BA.display, L = p === "none" ? "" : p)), BA.display = "inline-block")), D.overflow && (BA.overflow = "hidden", H.always(function() {
          BA.overflow = D.overflow[0], BA.overflowX = D.overflow[1], BA.overflowY = D.overflow[2];
        })), n = !1;
        for (w in u)
          n || (QA ? "hidden" in QA && (cA = QA.hidden) : QA = DA.access(Q, "fxshow", { display: L }), c && (QA.hidden = !cA), cA && Tg([Q], !0), H.done(function() {
            cA || Tg([Q]), DA.remove(Q, "fxshow");
            for (w in u)
              a.style(Q, w, u[w]);
          })), n = mB(cA ? QA[w] : 0, w, H), w in QA || (QA[w] = n.start, cA && (n.end = n.start, n.start = 0));
      }
    }
    function $Q(Q, i) {
      var D, w, e, c, G;
      for (D in Q)
        if (w = ag(D), e = i[w], c = Q[D], Array.isArray(c) && (e = c[1], c = Q[D] = c[0]), D !== w && (Q[w] = c, delete Q[D]), G = a.cssHooks[w], G && "expand" in G) {
          c = G.expand(c), delete Q[w];
          for (D in c)
            D in Q || (Q[D] = c[D], i[D] = e);
        } else
          i[w] = e;
    }
    function yg(Q, i, D) {
      var w, e, c = 0, G = yg.prefilters.length, R = a.Deferred().always(function() {
        delete n.elem;
      }), n = function() {
        if (e)
          return !1;
        for (var V = _g || xB(), H = Math.max(0, L.startTime + L.duration - V), u = H / L.duration || 0, BA = 1 - u, cA = 0, QA = L.tweens.length; cA < QA; cA++)
          L.tweens[cA].run(BA);
        return R.notifyWith(Q, [L, BA, H]), BA < 1 && QA ? H : (QA || R.notifyWith(Q, [L, 1, 0]), R.resolveWith(Q, [L]), !1);
      }, L = R.promise({
        elem: Q,
        props: a.extend({}, i),
        opts: a.extend(!0, {
          specialEasing: {},
          easing: a.easing._default
        }, D),
        originalProperties: i,
        originalOptions: D,
        startTime: _g || xB(),
        duration: D.duration,
        tweens: [],
        createTween: function(V, H) {
          var u = a.Tween(
            Q,
            L.opts,
            V,
            H,
            L.opts.specialEasing[V] || L.opts.easing
          );
          return L.tweens.push(u), u;
        },
        stop: function(V) {
          var H = 0, u = V ? L.tweens.length : 0;
          if (e)
            return this;
          for (e = !0; H < u; H++)
            L.tweens[H].run(1);
          return V ? (R.notifyWith(Q, [L, 1, 0]), R.resolveWith(Q, [L, V])) : R.rejectWith(Q, [L, V]), this;
        }
      }), p = L.props;
      for ($Q(p, L.opts.specialEasing); c < G; c++)
        if (w = yg.prefilters[c].call(L, Q, p, L.opts), w)
          return l(w.stop) && (a._queueHooks(L.elem, L.opts.queue).stop = w.stop.bind(w)), w;
      return a.map(p, mB, L), l(L.opts.start) && L.opts.start.call(Q, L), L.progress(L.opts.progress).done(L.opts.done, L.opts.complete).fail(L.opts.fail).always(L.opts.always), a.fx.timer(
        a.extend(n, {
          elem: Q,
          anim: L,
          queue: L.opts.queue
        })
      ), L;
    }
    a.Animation = a.extend(yg, {
      tweeners: {
        "*": [function(Q, i) {
          var D = this.createTween(Q, i);
          return NB(D.elem, Q, QI.exec(i), D), D;
        }]
      },
      tweener: function(Q, i) {
        l(Q) ? (i = Q, Q = ["*"]) : Q = Q.match(Dg);
        for (var D, w = 0, e = Q.length; w < e; w++)
          D = Q[w], yg.tweeners[D] = yg.tweeners[D] || [], yg.tweeners[D].unshift(i);
      },
      prefilters: [zQ],
      prefilter: function(Q, i) {
        i ? yg.prefilters.unshift(Q) : yg.prefilters.push(Q);
      }
    }), a.speed = function(Q, i, D) {
      var w = Q && typeof Q == "object" ? a.extend({}, Q) : {
        complete: D || !D && i || l(Q) && Q,
        duration: Q,
        easing: D && i || i && !l(i) && i
      };
      return a.fx.off ? w.duration = 0 : typeof w.duration != "number" && (w.duration in a.fx.speeds ? w.duration = a.fx.speeds[w.duration] : w.duration = a.fx.speeds._default), (w.queue == null || w.queue === !0) && (w.queue = "fx"), w.old = w.complete, w.complete = function() {
        l(w.old) && w.old.call(this), w.queue && a.dequeue(this, w.queue);
      }, w;
    }, a.fn.extend({
      fadeTo: function(Q, i, D, w) {
        return this.filter(KI).css("opacity", 0).show().end().animate({ opacity: i }, Q, D, w);
      },
      animate: function(Q, i, D, w) {
        var e = a.isEmptyObject(Q), c = a.speed(i, D, w), G = function() {
          var R = yg(this, a.extend({}, Q), c);
          (e || DA.get(this, "finish")) && R.stop(!0);
        };
        return G.finish = G, e || c.queue === !1 ? this.each(G) : this.queue(c.queue, G);
      },
      stop: function(Q, i, D) {
        var w = function(e) {
          var c = e.stop;
          delete e.stop, c(D);
        };
        return typeof Q != "string" && (D = i, i = Q, Q = void 0), i && this.queue(Q || "fx", []), this.each(function() {
          var e = !0, c = Q != null && Q + "queueHooks", G = a.timers, R = DA.get(this);
          if (c)
            R[c] && R[c].stop && w(R[c]);
          else
            for (c in R)
              R[c] && R[c].stop && XQ.test(c) && w(R[c]);
          for (c = G.length; c--; )
            G[c].elem === this && (Q == null || G[c].queue === Q) && (G[c].anim.stop(D), e = !1, G.splice(c, 1));
          (e || !D) && a.dequeue(this, Q);
        });
      },
      finish: function(Q) {
        return Q !== !1 && (Q = Q || "fx"), this.each(function() {
          var i, D = DA.get(this), w = D[Q + "queue"], e = D[Q + "queueHooks"], c = a.timers, G = w ? w.length : 0;
          for (D.finish = !0, a.queue(this, Q, []), e && e.stop && e.stop.call(this, !0), i = c.length; i--; )
            c[i].elem === this && c[i].queue === Q && (c[i].anim.stop(!0), c.splice(i, 1));
          for (i = 0; i < G; i++)
            w[i] && w[i].finish && w[i].finish.call(this);
          delete D.finish;
        });
      }
    }), a.each(["toggle", "show", "hide"], function(Q, i) {
      var D = a.fn[i];
      a.fn[i] = function(w, e, c) {
        return w == null || typeof w == "boolean" ? D.apply(this, arguments) : this.animate(UI(i, !0), w, e, c);
      };
    }), a.each({
      slideDown: UI("show"),
      slideUp: UI("hide"),
      slideToggle: UI("toggle"),
      fadeIn: { opacity: "show" },
      fadeOut: { opacity: "hide" },
      fadeToggle: { opacity: "toggle" }
    }, function(Q, i) {
      a.fn[Q] = function(D, w, e) {
        return this.animate(i, D, w, e);
      };
    }), a.timers = [], a.fx.tick = function() {
      var Q, i = 0, D = a.timers;
      for (_g = Date.now(); i < D.length; i++)
        Q = D[i], !Q() && D[i] === Q && D.splice(i--, 1);
      D.length || a.fx.stop(), _g = void 0;
    }, a.fx.timer = function(Q) {
      a.timers.push(Q), a.fx.start();
    }, a.fx.interval = 13, a.fx.start = function() {
      JI || (JI = !0, PI());
    }, a.fx.stop = function() {
      JI = null;
    }, a.fx.speeds = {
      slow: 600,
      fast: 200,
      // Default speed
      _default: 400
    }, a.fn.delay = function(Q, i) {
      return Q = a.fx && a.fx.speeds[Q] || Q, i = i || "fx", this.queue(i, function(D, w) {
        var e = A.setTimeout(D, Q);
        w.stop = function() {
          A.clearTimeout(e);
        };
      });
    }, function() {
      var Q = q.createElement("input"), i = q.createElement("select"), D = i.appendChild(q.createElement("option"));
      Q.type = "checkbox", K.checkOn = Q.value !== "", K.optSelected = D.selected, Q = q.createElement("input"), Q.value = "t", Q.type = "radio", K.radioValue = Q.value === "t";
    }();
    var OB, iI = a.expr.attrHandle;
    a.fn.extend({
      attr: function(Q, i) {
        return HA(this, a.attr, Q, i, arguments.length > 1);
      },
      removeAttr: function(Q) {
        return this.each(function() {
          a.removeAttr(this, Q);
        });
      }
    }), a.extend({
      attr: function(Q, i, D) {
        var w, e, c = Q.nodeType;
        if (!(c === 3 || c === 8 || c === 2)) {
          if (typeof Q.getAttribute > "u")
            return a.prop(Q, i, D);
          if ((c !== 1 || !a.isXMLDoc(Q)) && (e = a.attrHooks[i.toLowerCase()] || (a.expr.match.bool.test(i) ? OB : void 0)), D !== void 0) {
            if (D === null) {
              a.removeAttr(Q, i);
              return;
            }
            return e && "set" in e && (w = e.set(Q, D, i)) !== void 0 ? w : (Q.setAttribute(i, D + ""), D);
          }
          return e && "get" in e && (w = e.get(Q, i)) !== null ? w : (w = a.find.attr(Q, i), w ?? void 0);
        }
      },
      attrHooks: {
        type: {
          set: function(Q, i) {
            if (!K.radioValue && i === "radio" && yA(Q, "input")) {
              var D = Q.value;
              return Q.setAttribute("type", i), D && (Q.value = D), i;
            }
          }
        }
      },
      removeAttr: function(Q, i) {
        var D, w = 0, e = i && i.match(Dg);
        if (e && Q.nodeType === 1)
          for (; D = e[w++]; )
            Q.removeAttribute(D);
      }
    }), OB = {
      set: function(Q, i, D) {
        return i === !1 ? a.removeAttr(Q, D) : Q.setAttribute(D, D), D;
      }
    }, a.each(a.expr.match.bool.source.match(/\w+/g), function(Q, i) {
      var D = iI[i] || a.find.attr;
      iI[i] = function(w, e, c) {
        var G, R, n = e.toLowerCase();
        return c || (R = iI[n], iI[n] = G, G = D(w, e, c) != null ? n : null, iI[n] = R), G;
      };
    });
    var AC = /^(?:input|select|textarea|button)$/i, gC = /^(?:a|area)$/i;
    a.fn.extend({
      prop: function(Q, i) {
        return HA(this, a.prop, Q, i, arguments.length > 1);
      },
      removeProp: function(Q) {
        return this.each(function() {
          delete this[a.propFix[Q] || Q];
        });
      }
    }), a.extend({
      prop: function(Q, i, D) {
        var w, e, c = Q.nodeType;
        if (!(c === 3 || c === 8 || c === 2))
          return (c !== 1 || !a.isXMLDoc(Q)) && (i = a.propFix[i] || i, e = a.propHooks[i]), D !== void 0 ? e && "set" in e && (w = e.set(Q, D, i)) !== void 0 ? w : Q[i] = D : e && "get" in e && (w = e.get(Q, i)) !== null ? w : Q[i];
      },
      propHooks: {
        tabIndex: {
          get: function(Q) {
            var i = a.find.attr(Q, "tabindex");
            return i ? parseInt(i, 10) : AC.test(Q.nodeName) || gC.test(Q.nodeName) && Q.href ? 0 : -1;
          }
        }
      },
      propFix: {
        for: "htmlFor",
        class: "className"
      }
    }), K.optSelected || (a.propHooks.selected = {
      get: function(Q) {
        var i = Q.parentNode;
        return i && i.parentNode && i.parentNode.selectedIndex, null;
      },
      set: function(Q) {
        var i = Q.parentNode;
        i && (i.selectedIndex, i.parentNode && i.parentNode.selectedIndex);
      }
    }), a.each([
      "tabIndex",
      "readOnly",
      "maxLength",
      "cellSpacing",
      "cellPadding",
      "rowSpan",
      "colSpan",
      "useMap",
      "frameBorder",
      "contentEditable"
    ], function() {
      a.propFix[this.toLowerCase()] = this;
    });
    function Hg(Q) {
      var i = Q.match(Dg) || [];
      return i.join(" ");
    }
    function pg(Q) {
      return Q.getAttribute && Q.getAttribute("class") || "";
    }
    function XI(Q) {
      return Array.isArray(Q) ? Q : typeof Q == "string" ? Q.match(Dg) || [] : [];
    }
    a.fn.extend({
      addClass: function(Q) {
        var i, D, w, e, c, G;
        return l(Q) ? this.each(function(R) {
          a(this).addClass(Q.call(this, R, pg(this)));
        }) : (i = XI(Q), i.length ? this.each(function() {
          if (w = pg(this), D = this.nodeType === 1 && " " + Hg(w) + " ", D) {
            for (c = 0; c < i.length; c++)
              e = i[c], D.indexOf(" " + e + " ") < 0 && (D += e + " ");
            G = Hg(D), w !== G && this.setAttribute("class", G);
          }
        }) : this);
      },
      removeClass: function(Q) {
        var i, D, w, e, c, G;
        return l(Q) ? this.each(function(R) {
          a(this).removeClass(Q.call(this, R, pg(this)));
        }) : arguments.length ? (i = XI(Q), i.length ? this.each(function() {
          if (w = pg(this), D = this.nodeType === 1 && " " + Hg(w) + " ", D) {
            for (c = 0; c < i.length; c++)
              for (e = i[c]; D.indexOf(" " + e + " ") > -1; )
                D = D.replace(" " + e + " ", " ");
            G = Hg(D), w !== G && this.setAttribute("class", G);
          }
        }) : this) : this.attr("class", "");
      },
      toggleClass: function(Q, i) {
        var D, w, e, c, G = typeof Q, R = G === "string" || Array.isArray(Q);
        return l(Q) ? this.each(function(n) {
          a(this).toggleClass(
            Q.call(this, n, pg(this), i),
            i
          );
        }) : typeof i == "boolean" && R ? i ? this.addClass(Q) : this.removeClass(Q) : (D = XI(Q), this.each(function() {
          if (R)
            for (c = a(this), e = 0; e < D.length; e++)
              w = D[e], c.hasClass(w) ? c.removeClass(w) : c.addClass(w);
          else
            (Q === void 0 || G === "boolean") && (w = pg(this), w && DA.set(this, "__className__", w), this.setAttribute && this.setAttribute(
              "class",
              w || Q === !1 ? "" : DA.get(this, "__className__") || ""
            ));
        }));
      },
      hasClass: function(Q) {
        var i, D, w = 0;
        for (i = " " + Q + " "; D = this[w++]; )
          if (D.nodeType === 1 && (" " + Hg(pg(D)) + " ").indexOf(i) > -1)
            return !0;
        return !1;
      }
    });
    var IC = /\r/g;
    a.fn.extend({
      val: function(Q) {
        var i, D, w, e = this[0];
        return arguments.length ? (w = l(Q), this.each(function(c) {
          var G;
          this.nodeType === 1 && (w ? G = Q.call(this, c, a(this).val()) : G = Q, G == null ? G = "" : typeof G == "number" ? G += "" : Array.isArray(G) && (G = a.map(G, function(R) {
            return R == null ? "" : R + "";
          })), i = a.valHooks[this.type] || a.valHooks[this.nodeName.toLowerCase()], (!i || !("set" in i) || i.set(this, G, "value") === void 0) && (this.value = G));
        })) : e ? (i = a.valHooks[e.type] || a.valHooks[e.nodeName.toLowerCase()], i && "get" in i && (D = i.get(e, "value")) !== void 0 ? D : (D = e.value, typeof D == "string" ? D.replace(IC, "") : D ?? "")) : void 0;
      }
    }), a.extend({
      valHooks: {
        option: {
          get: function(Q) {
            var i = a.find.attr(Q, "value");
            return i ?? // Support: IE <=10 - 11 only
            // option.text throws exceptions (trac-14686, trac-14858)
            // Strip and collapse whitespace
            // https://html.spec.whatwg.org/#strip-and-collapse-whitespace
            Hg(a.text(Q));
          }
        },
        select: {
          get: function(Q) {
            var i, D, w, e = Q.options, c = Q.selectedIndex, G = Q.type === "select-one", R = G ? null : [], n = G ? c + 1 : e.length;
            for (c < 0 ? w = n : w = G ? c : 0; w < n; w++)
              if (D = e[w], (D.selected || w === c) && // Don't return options that are disabled or in a disabled optgroup
              !D.disabled && (!D.parentNode.disabled || !yA(D.parentNode, "optgroup"))) {
                if (i = a(D).val(), G)
                  return i;
                R.push(i);
              }
            return R;
          },
          set: function(Q, i) {
            for (var D, w, e = Q.options, c = a.makeArray(i), G = e.length; G--; )
              w = e[G], (w.selected = a.inArray(a.valHooks.option.get(w), c) > -1) && (D = !0);
            return D || (Q.selectedIndex = -1), c;
          }
        }
      }
    }), a.each(["radio", "checkbox"], function() {
      a.valHooks[this] = {
        set: function(Q, i) {
          if (Array.isArray(i))
            return Q.checked = a.inArray(a(Q).val(), i) > -1;
        }
      }, K.checkOn || (a.valHooks[this].get = function(Q) {
        return Q.getAttribute("value") === null ? "on" : Q.value;
      });
    }), K.focusin = "onfocusin" in A;
    var vB = /^(?:focusinfocus|focusoutblur)$/, bB = function(Q) {
      Q.stopPropagation();
    };
    a.extend(a.event, {
      trigger: function(Q, i, D, w) {
        var e, c, G, R, n, L, p, V, H = [D || q], u = k.call(Q, "type") ? Q.type : Q, BA = k.call(Q, "namespace") ? Q.namespace.split(".") : [];
        if (c = V = G = D = D || q, !(D.nodeType === 3 || D.nodeType === 8) && !vB.test(u + a.event.triggered) && (u.indexOf(".") > -1 && (BA = u.split("."), u = BA.shift(), BA.sort()), n = u.indexOf(":") < 0 && "on" + u, Q = Q[a.expando] ? Q : new a.Event(u, typeof Q == "object" && Q), Q.isTrigger = w ? 2 : 3, Q.namespace = BA.join("."), Q.rnamespace = Q.namespace ? new RegExp("(^|\\.)" + BA.join("\\.(?:.*\\.|)") + "(\\.|$)") : null, Q.result = void 0, Q.target || (Q.target = D), i = i == null ? [Q] : a.makeArray(i, [Q]), p = a.event.special[u] || {}, !(!w && p.trigger && p.trigger.apply(D, i) === !1))) {
          if (!w && !p.noBubble && !d(D)) {
            for (R = p.delegateType || u, vB.test(R + u) || (c = c.parentNode); c; c = c.parentNode)
              H.push(c), G = c;
            G === (D.ownerDocument || q) && H.push(G.defaultView || G.parentWindow || A);
          }
          for (e = 0; (c = H[e++]) && !Q.isPropagationStopped(); )
            V = c, Q.type = e > 1 ? R : p.bindType || u, L = (DA.get(c, "events") || /* @__PURE__ */ Object.create(null))[Q.type] && DA.get(c, "handle"), L && L.apply(c, i), L = n && c[n], L && L.apply && II(c) && (Q.result = L.apply(c, i), Q.result === !1 && Q.preventDefault());
          return Q.type = u, !w && !Q.isDefaultPrevented() && (!p._default || p._default.apply(H.pop(), i) === !1) && II(D) && n && l(D[u]) && !d(D) && (G = D[n], G && (D[n] = null), a.event.triggered = u, Q.isPropagationStopped() && V.addEventListener(u, bB), D[u](), Q.isPropagationStopped() && V.removeEventListener(u, bB), a.event.triggered = void 0, G && (D[n] = G)), Q.result;
        }
      },
      // Piggyback on a donor event to simulate a different one
      // Used only for `focus(in | out)` events
      simulate: function(Q, i, D) {
        var w = a.extend(
          new a.Event(),
          D,
          {
            type: Q,
            isSimulated: !0
          }
        );
        a.event.trigger(w, null, i);
      }
    }), a.fn.extend({
      trigger: function(Q, i) {
        return this.each(function() {
          a.event.trigger(Q, i, this);
        });
      },
      triggerHandler: function(Q, i) {
        var D = this[0];
        if (D)
          return a.event.trigger(Q, i, D, !0);
      }
    }), K.focusin || a.each({ focus: "focusin", blur: "focusout" }, function(Q, i) {
      var D = function(w) {
        a.event.simulate(i, w.target, a.event.fix(w));
      };
      a.event.special[i] = {
        setup: function() {
          var w = this.ownerDocument || this.document || this, e = DA.access(w, i);
          e || w.addEventListener(Q, D, !0), DA.access(w, i, (e || 0) + 1);
        },
        teardown: function() {
          var w = this.ownerDocument || this.document || this, e = DA.access(w, i) - 1;
          e ? DA.access(w, i, e) : (w.removeEventListener(Q, D, !0), DA.remove(w, i));
        }
      };
    });
    var oI = A.location, ZB = { guid: Date.now() }, zI = /\?/;
    a.parseXML = function(Q) {
      var i, D;
      if (!Q || typeof Q != "string")
        return null;
      try {
        i = new A.DOMParser().parseFromString(Q, "text/xml");
      } catch {
      }
      return D = i && i.getElementsByTagName("parsererror")[0], (!i || D) && a.error("Invalid XML: " + (D ? a.map(D.childNodes, function(w) {
        return w.textContent;
      }).join(`
`) : Q)), i;
    };
    var BC = /\[\]$/, TB = /\r?\n/g, QC = /^(?:submit|button|image|reset|file)$/i, CC = /^(?:input|select|textarea|keygen)/i;
    function $I(Q, i, D, w) {
      var e;
      if (Array.isArray(i))
        a.each(i, function(c, G) {
          D || BC.test(Q) ? w(Q, G) : $I(
            Q + "[" + (typeof G == "object" && G != null ? c : "") + "]",
            G,
            D,
            w
          );
        });
      else if (!D && v(i) === "object")
        for (e in i)
          $I(Q + "[" + e + "]", i[e], D, w);
      else
        w(Q, i);
    }
    a.param = function(Q, i) {
      var D, w = [], e = function(c, G) {
        var R = l(G) ? G() : G;
        w[w.length] = encodeURIComponent(c) + "=" + encodeURIComponent(R ?? "");
      };
      if (Q == null)
        return "";
      if (Array.isArray(Q) || Q.jquery && !a.isPlainObject(Q))
        a.each(Q, function() {
          e(this.name, this.value);
        });
      else
        for (D in Q)
          $I(D, Q[D], i, e);
      return w.join("&");
    }, a.fn.extend({
      serialize: function() {
        return a.param(this.serializeArray());
      },
      serializeArray: function() {
        return this.map(function() {
          var Q = a.prop(this, "elements");
          return Q ? a.makeArray(Q) : this;
        }).filter(function() {
          var Q = this.type;
          return this.name && !a(this).is(":disabled") && CC.test(this.nodeName) && !QC.test(Q) && (this.checked || !CI.test(Q));
        }).map(function(Q, i) {
          var D = a(this).val();
          return D == null ? null : Array.isArray(D) ? a.map(D, function(w) {
            return { name: i.name, value: w.replace(TB, `\r
`) };
          }) : { name: i.name, value: D.replace(TB, `\r
`) };
        }).get();
      }
    });
    var EC = /%20/g, iC = /#.*$/, oC = /([?&])_=[^&]*/, DC = /^(.*?):[ \t]*([^\r\n]*)$/mg, aC = /^(?:about|app|app-storage|.+-extension|file|res|widget):$/, wC = /^(?:GET|HEAD)$/, tC = /^\/\//, WB = {}, AB = {}, VB = "*/".concat("*"), gB = q.createElement("a");
    gB.href = oI.href;
    function jB(Q) {
      return function(i, D) {
        typeof i != "string" && (D = i, i = "*");
        var w, e = 0, c = i.toLowerCase().match(Dg) || [];
        if (l(D))
          for (; w = c[e++]; )
            w[0] === "+" ? (w = w.slice(1) || "*", (Q[w] = Q[w] || []).unshift(D)) : (Q[w] = Q[w] || []).push(D);
      };
    }
    function _B(Q, i, D, w) {
      var e = {}, c = Q === AB;
      function G(R) {
        var n;
        return e[R] = !0, a.each(Q[R] || [], function(L, p) {
          var V = p(i, D, w);
          if (typeof V == "string" && !c && !e[V])
            return i.dataTypes.unshift(V), G(V), !1;
          if (c)
            return !(n = V);
        }), n;
      }
      return G(i.dataTypes[0]) || !e["*"] && G("*");
    }
    function IB(Q, i) {
      var D, w, e = a.ajaxSettings.flatOptions || {};
      for (D in i)
        i[D] !== void 0 && ((e[D] ? Q : w || (w = {}))[D] = i[D]);
      return w && a.extend(!0, Q, w), Q;
    }
    function sC(Q, i, D) {
      for (var w, e, c, G, R = Q.contents, n = Q.dataTypes; n[0] === "*"; )
        n.shift(), w === void 0 && (w = Q.mimeType || i.getResponseHeader("Content-Type"));
      if (w) {
        for (e in R)
          if (R[e] && R[e].test(w)) {
            n.unshift(e);
            break;
          }
      }
      if (n[0] in D)
        c = n[0];
      else {
        for (e in D) {
          if (!n[0] || Q.converters[e + " " + n[0]]) {
            c = e;
            break;
          }
          G || (G = e);
        }
        c = c || G;
      }
      if (c)
        return c !== n[0] && n.unshift(c), D[c];
    }
    function eC(Q, i, D, w) {
      var e, c, G, R, n, L = {}, p = Q.dataTypes.slice();
      if (p[1])
        for (G in Q.converters)
          L[G.toLowerCase()] = Q.converters[G];
      for (c = p.shift(); c; )
        if (Q.responseFields[c] && (D[Q.responseFields[c]] = i), !n && w && Q.dataFilter && (i = Q.dataFilter(i, Q.dataType)), n = c, c = p.shift(), c) {
          if (c === "*")
            c = n;
          else if (n !== "*" && n !== c) {
            if (G = L[n + " " + c] || L["* " + c], !G) {
              for (e in L)
                if (R = e.split(" "), R[1] === c && (G = L[n + " " + R[0]] || L["* " + R[0]], G)) {
                  G === !0 ? G = L[e] : L[e] !== !0 && (c = R[0], p.unshift(R[1]));
                  break;
                }
            }
            if (G !== !0)
              if (G && Q.throws)
                i = G(i);
              else
                try {
                  i = G(i);
                } catch (V) {
                  return {
                    state: "parsererror",
                    error: G ? V : "No conversion from " + n + " to " + c
                  };
                }
          }
        }
      return { state: "success", data: i };
    }
    a.extend({
      // Counter for holding the number of active queries
      active: 0,
      // Last-Modified header cache for next request
      lastModified: {},
      etag: {},
      ajaxSettings: {
        url: oI.href,
        type: "GET",
        isLocal: aC.test(oI.protocol),
        global: !0,
        processData: !0,
        async: !0,
        contentType: "application/x-www-form-urlencoded; charset=UTF-8",
        /*
        timeout: 0,
        data: null,
        dataType: null,
        username: null,
        password: null,
        cache: null,
        throws: false,
        traditional: false,
        headers: {},
        */
        accepts: {
          "*": VB,
          text: "text/plain",
          html: "text/html",
          xml: "application/xml, text/xml",
          json: "application/json, text/javascript"
        },
        contents: {
          xml: /\bxml\b/,
          html: /\bhtml/,
          json: /\bjson\b/
        },
        responseFields: {
          xml: "responseXML",
          text: "responseText",
          json: "responseJSON"
        },
        // Data converters
        // Keys separate source (or catchall "*") and destination types with a single space
        converters: {
          // Convert anything to text
          "* text": String,
          // Text to html (true = no transformation)
          "text html": !0,
          // Evaluate text as a json expression
          "text json": JSON.parse,
          // Parse text as xml
          "text xml": a.parseXML
        },
        // For options that shouldn't be deep extended:
        // you can add your own custom options here if
        // and when you create one that shouldn't be
        // deep extended (see ajaxExtend)
        flatOptions: {
          url: !0,
          context: !0
        }
      },
      // Creates a full fledged settings object into target
      // with both ajaxSettings and settings fields.
      // If target is omitted, writes into ajaxSettings.
      ajaxSetup: function(Q, i) {
        return i ? (
          // Building a settings object
          IB(IB(Q, a.ajaxSettings), i)
        ) : (
          // Extending ajaxSettings
          IB(a.ajaxSettings, Q)
        );
      },
      ajaxPrefilter: jB(WB),
      ajaxTransport: jB(AB),
      // Main method
      ajax: function(Q, i) {
        typeof Q == "object" && (i = Q, Q = void 0), i = i || {};
        var D, w, e, c, G, R, n, L, p, V, H = a.ajaxSetup({}, i), u = H.context || H, BA = H.context && (u.nodeType || u.jquery) ? a(u) : a.event, cA = a.Deferred(), QA = a.Callbacks("once memory"), TA = H.statusCode || {}, bA = {}, wg = {}, UA = "canceled", eA = {
          readyState: 0,
          // Builds headers hashtable if needed
          getResponseHeader: function(YA) {
            var fA;
            if (n) {
              if (!c)
                for (c = {}; fA = DC.exec(e); )
                  c[fA[1].toLowerCase() + " "] = (c[fA[1].toLowerCase() + " "] || []).concat(fA[2]);
              fA = c[YA.toLowerCase() + " "];
            }
            return fA == null ? null : fA.join(", ");
          },
          // Raw string
          getAllResponseHeaders: function() {
            return n ? e : null;
          },
          // Caches the header
          setRequestHeader: function(YA, fA) {
            return n == null && (YA = wg[YA.toLowerCase()] = wg[YA.toLowerCase()] || YA, bA[YA] = fA), this;
          },
          // Overrides response content-type header
          overrideMimeType: function(YA) {
            return n == null && (H.mimeType = YA), this;
          },
          // Status-dependent callbacks
          statusCode: function(YA) {
            var fA;
            if (YA)
              if (n)
                eA.always(YA[eA.status]);
              else
                for (fA in YA)
                  TA[fA] = [TA[fA], YA[fA]];
            return this;
          },
          // Cancel the request
          abort: function(YA) {
            var fA = YA || UA;
            return D && D.abort(fA), Eg(0, fA), this;
          }
        };
        if (cA.promise(eA), H.url = ((Q || H.url || oI.href) + "").replace(tC, oI.protocol + "//"), H.type = i.method || i.type || H.method || H.type, H.dataTypes = (H.dataType || "*").toLowerCase().match(Dg) || [""], H.crossDomain == null) {
          R = q.createElement("a");
          try {
            R.href = H.url, R.href = R.href, H.crossDomain = gB.protocol + "//" + gB.host != R.protocol + "//" + R.host;
          } catch {
            H.crossDomain = !0;
          }
        }
        if (H.data && H.processData && typeof H.data != "string" && (H.data = a.param(H.data, H.traditional)), _B(WB, H, i, eA), n)
          return eA;
        L = a.event && H.global, L && a.active++ === 0 && a.event.trigger("ajaxStart"), H.type = H.type.toUpperCase(), H.hasContent = !wC.test(H.type), w = H.url.replace(iC, ""), H.hasContent ? H.data && H.processData && (H.contentType || "").indexOf("application/x-www-form-urlencoded") === 0 && (H.data = H.data.replace(EC, "+")) : (V = H.url.slice(w.length), H.data && (H.processData || typeof H.data == "string") && (w += (zI.test(w) ? "&" : "?") + H.data, delete H.data), H.cache === !1 && (w = w.replace(oC, "$1"), V = (zI.test(w) ? "&" : "?") + "_=" + ZB.guid++ + V), H.url = w + V), H.ifModified && (a.lastModified[w] && eA.setRequestHeader("If-Modified-Since", a.lastModified[w]), a.etag[w] && eA.setRequestHeader("If-None-Match", a.etag[w])), (H.data && H.hasContent && H.contentType !== !1 || i.contentType) && eA.setRequestHeader("Content-Type", H.contentType), eA.setRequestHeader(
          "Accept",
          H.dataTypes[0] && H.accepts[H.dataTypes[0]] ? H.accepts[H.dataTypes[0]] + (H.dataTypes[0] !== "*" ? ", " + VB + "; q=0.01" : "") : H.accepts["*"]
        );
        for (p in H.headers)
          eA.setRequestHeader(p, H.headers[p]);
        if (H.beforeSend && (H.beforeSend.call(u, eA, H) === !1 || n))
          return eA.abort();
        if (UA = "abort", QA.add(H.complete), eA.done(H.success), eA.fail(H.error), D = _B(AB, H, i, eA), !D)
          Eg(-1, "No Transport");
        else {
          if (eA.readyState = 1, L && BA.trigger("ajaxSend", [eA, H]), n)
            return eA;
          H.async && H.timeout > 0 && (G = A.setTimeout(function() {
            eA.abort("timeout");
          }, H.timeout));
          try {
            n = !1, D.send(bA, Eg);
          } catch (YA) {
            if (n)
              throw YA;
            Eg(-1, YA);
          }
        }
        function Eg(YA, fA, aI, SI) {
          var tg, fg, ug, ig, Ug, hg = fA;
          n || (n = !0, G && A.clearTimeout(G), D = void 0, e = SI || "", eA.readyState = YA > 0 ? 4 : 0, tg = YA >= 200 && YA < 300 || YA === 304, aI && (ig = sC(H, eA, aI)), !tg && a.inArray("script", H.dataTypes) > -1 && a.inArray("json", H.dataTypes) < 0 && (H.converters["text script"] = function() {
          }), ig = eC(H, ig, eA, tg), tg ? (H.ifModified && (Ug = eA.getResponseHeader("Last-Modified"), Ug && (a.lastModified[w] = Ug), Ug = eA.getResponseHeader("etag"), Ug && (a.etag[w] = Ug)), YA === 204 || H.type === "HEAD" ? hg = "nocontent" : YA === 304 ? hg = "notmodified" : (hg = ig.state, fg = ig.data, ug = ig.error, tg = !ug)) : (ug = hg, (YA || !hg) && (hg = "error", YA < 0 && (YA = 0))), eA.status = YA, eA.statusText = (fA || hg) + "", tg ? cA.resolveWith(u, [fg, hg, eA]) : cA.rejectWith(u, [eA, hg, ug]), eA.statusCode(TA), TA = void 0, L && BA.trigger(
            tg ? "ajaxSuccess" : "ajaxError",
            [eA, H, tg ? fg : ug]
          ), QA.fireWith(u, [eA, hg]), L && (BA.trigger("ajaxComplete", [eA, H]), --a.active || a.event.trigger("ajaxStop")));
        }
        return eA;
      },
      getJSON: function(Q, i, D) {
        return a.get(Q, i, D, "json");
      },
      getScript: function(Q, i) {
        return a.get(Q, void 0, i, "script");
      }
    }), a.each(["get", "post"], function(Q, i) {
      a[i] = function(D, w, e, c) {
        return l(w) && (c = c || e, e = w, w = void 0), a.ajax(a.extend({
          url: D,
          type: i,
          dataType: c,
          data: w,
          success: e
        }, a.isPlainObject(D) && D));
      };
    }), a.ajaxPrefilter(function(Q) {
      var i;
      for (i in Q.headers)
        i.toLowerCase() === "content-type" && (Q.contentType = Q.headers[i] || "");
    }), a._evalUrl = function(Q, i, D) {
      return a.ajax({
        url: Q,
        // Make this explicit, since user can override this through ajaxSetup (trac-11264)
        type: "GET",
        dataType: "script",
        cache: !0,
        async: !1,
        global: !1,
        // Only evaluate the response if it is successful (gh-4126)
        // dataFilter is not invoked for failure responses, so using it instead
        // of the default converter is kludgy but it works.
        converters: {
          "text script": function() {
          }
        },
        dataFilter: function(w) {
          a.globalEval(w, i, D);
        }
      });
    }, a.fn.extend({
      wrapAll: function(Q) {
        var i;
        return this[0] && (l(Q) && (Q = Q.call(this[0])), i = a(Q, this[0].ownerDocument).eq(0).clone(!0), this[0].parentNode && i.insertBefore(this[0]), i.map(function() {
          for (var D = this; D.firstElementChild; )
            D = D.firstElementChild;
          return D;
        }).append(this)), this;
      },
      wrapInner: function(Q) {
        return l(Q) ? this.each(function(i) {
          a(this).wrapInner(Q.call(this, i));
        }) : this.each(function() {
          var i = a(this), D = i.contents();
          D.length ? D.wrapAll(Q) : i.append(Q);
        });
      },
      wrap: function(Q) {
        var i = l(Q);
        return this.each(function(D) {
          a(this).wrapAll(i ? Q.call(this, D) : Q);
        });
      },
      unwrap: function(Q) {
        return this.parent(Q).not("body").each(function() {
          a(this).replaceWith(this.childNodes);
        }), this;
      }
    }), a.expr.pseudos.hidden = function(Q) {
      return !a.expr.pseudos.visible(Q);
    }, a.expr.pseudos.visible = function(Q) {
      return !!(Q.offsetWidth || Q.offsetHeight || Q.getClientRects().length);
    }, a.ajaxSettings.xhr = function() {
      try {
        return new A.XMLHttpRequest();
      } catch {
      }
    };
    var cC = {
      // File protocol always yields status code 0, assume 200
      0: 200,
      // Support: IE <=9 only
      // trac-1450: sometimes IE returns 1223 when it should be 204
      1223: 204
    }, DI = a.ajaxSettings.xhr();
    K.cors = !!DI && "withCredentials" in DI, K.ajax = DI = !!DI, a.ajaxTransport(function(Q) {
      var i, D;
      if (K.cors || DI && !Q.crossDomain)
        return {
          send: function(w, e) {
            var c, G = Q.xhr();
            if (G.open(
              Q.type,
              Q.url,
              Q.async,
              Q.username,
              Q.password
            ), Q.xhrFields)
              for (c in Q.xhrFields)
                G[c] = Q.xhrFields[c];
            Q.mimeType && G.overrideMimeType && G.overrideMimeType(Q.mimeType), !Q.crossDomain && !w["X-Requested-With"] && (w["X-Requested-With"] = "XMLHttpRequest");
            for (c in w)
              G.setRequestHeader(c, w[c]);
            i = function(R) {
              return function() {
                i && (i = D = G.onload = G.onerror = G.onabort = G.ontimeout = G.onreadystatechange = null, R === "abort" ? G.abort() : R === "error" ? typeof G.status != "number" ? e(0, "error") : e(
                  // File: protocol always yields status 0; see trac-8605, trac-14207
                  G.status,
                  G.statusText
                ) : e(
                  cC[G.status] || G.status,
                  G.statusText,
                  // Support: IE <=9 only
                  // IE9 has no XHR2 but throws on binary (trac-11426)
                  // For XHR2 non-text, let the caller handle it (gh-2498)
                  (G.responseType || "text") !== "text" || typeof G.responseText != "string" ? { binary: G.response } : { text: G.responseText },
                  G.getAllResponseHeaders()
                ));
              };
            }, G.onload = i(), D = G.onerror = G.ontimeout = i("error"), G.onabort !== void 0 ? G.onabort = D : G.onreadystatechange = function() {
              G.readyState === 4 && A.setTimeout(function() {
                i && D();
              });
            }, i = i("abort");
            try {
              G.send(Q.hasContent && Q.data || null);
            } catch (R) {
              if (i)
                throw R;
            }
          },
          abort: function() {
            i && i();
          }
        };
    }), a.ajaxPrefilter(function(Q) {
      Q.crossDomain && (Q.contents.script = !1);
    }), a.ajaxSetup({
      accepts: {
        script: "text/javascript, application/javascript, application/ecmascript, application/x-ecmascript"
      },
      contents: {
        script: /\b(?:java|ecma)script\b/
      },
      converters: {
        "text script": function(Q) {
          return a.globalEval(Q), Q;
        }
      }
    }), a.ajaxPrefilter("script", function(Q) {
      Q.cache === void 0 && (Q.cache = !1), Q.crossDomain && (Q.type = "GET");
    }), a.ajaxTransport("script", function(Q) {
      if (Q.crossDomain || Q.scriptAttrs) {
        var i, D;
        return {
          send: function(w, e) {
            i = a("<script>").attr(Q.scriptAttrs || {}).prop({ charset: Q.scriptCharset, src: Q.url }).on("load error", D = function(c) {
              i.remove(), D = null, c && e(c.type === "error" ? 404 : 200, c.type);
            }), q.head.appendChild(i[0]);
          },
          abort: function() {
            D && D();
          }
        };
      }
    });
    var PB = [], BB = /(=)\?(?=&|$)|\?\?/;
    a.ajaxSetup({
      jsonp: "callback",
      jsonpCallback: function() {
        var Q = PB.pop() || a.expando + "_" + ZB.guid++;
        return this[Q] = !0, Q;
      }
    }), a.ajaxPrefilter("json jsonp", function(Q, i, D) {
      var w, e, c, G = Q.jsonp !== !1 && (BB.test(Q.url) ? "url" : typeof Q.data == "string" && (Q.contentType || "").indexOf("application/x-www-form-urlencoded") === 0 && BB.test(Q.data) && "data");
      if (G || Q.dataTypes[0] === "jsonp")
        return w = Q.jsonpCallback = l(Q.jsonpCallback) ? Q.jsonpCallback() : Q.jsonpCallback, G ? Q[G] = Q[G].replace(BB, "$1" + w) : Q.jsonp !== !1 && (Q.url += (zI.test(Q.url) ? "&" : "?") + Q.jsonp + "=" + w), Q.converters["script json"] = function() {
          return c || a.error(w + " was not called"), c[0];
        }, Q.dataTypes[0] = "json", e = A[w], A[w] = function() {
          c = arguments;
        }, D.always(function() {
          e === void 0 ? a(A).removeProp(w) : A[w] = e, Q[w] && (Q.jsonpCallback = i.jsonpCallback, PB.push(w)), c && l(e) && e(c[0]), c = e = void 0;
        }), "script";
    }), K.createHTMLDocument = function() {
      var Q = q.implementation.createHTMLDocument("").body;
      return Q.innerHTML = "<form></form><form></form>", Q.childNodes.length === 2;
    }(), a.parseHTML = function(Q, i, D) {
      if (typeof Q != "string")
        return [];
      typeof i == "boolean" && (D = i, i = !1);
      var w, e, c;
      return i || (K.createHTMLDocument ? (i = q.implementation.createHTMLDocument(""), w = i.createElement("base"), w.href = q.location.href, i.head.appendChild(w)) : i = q), e = rA.exec(Q), c = !D && [], e ? [i.createElement(e[1])] : (e = RB([Q], i, c), c && c.length && a(c).remove(), a.merge([], e.childNodes));
    }, a.fn.load = function(Q, i, D) {
      var w, e, c, G = this, R = Q.indexOf(" ");
      return R > -1 && (w = Hg(Q.slice(R)), Q = Q.slice(0, R)), l(i) ? (D = i, i = void 0) : i && typeof i == "object" && (e = "POST"), G.length > 0 && a.ajax({
        url: Q,
        // If "type" variable is undefined, then "GET" method will be used.
        // Make value of this field explicit since
        // user can override it through ajaxSetup method
        type: e || "GET",
        dataType: "html",
        data: i
      }).done(function(n) {
        c = arguments, G.html(w ? (
          // If a selector was specified, locate the right elements in a dummy div
          // Exclude scripts to avoid IE 'Permission Denied' errors
          a("<div>").append(a.parseHTML(n)).find(w)
        ) : (
          // Otherwise use the full result
          n
        ));
      }).always(D && function(n, L) {
        G.each(function() {
          D.apply(this, c || [n.responseText, L, n]);
        });
      }), this;
    }, a.expr.pseudos.animated = function(Q) {
      return a.grep(a.timers, function(i) {
        return Q === i.elem;
      }).length;
    }, a.offset = {
      setOffset: function(Q, i, D) {
        var w, e, c, G, R, n, L, p = a.css(Q, "position"), V = a(Q), H = {};
        p === "static" && (Q.style.position = "relative"), R = V.offset(), c = a.css(Q, "top"), n = a.css(Q, "left"), L = (p === "absolute" || p === "fixed") && (c + n).indexOf("auto") > -1, L ? (w = V.position(), G = w.top, e = w.left) : (G = parseFloat(c) || 0, e = parseFloat(n) || 0), l(i) && (i = i.call(Q, D, a.extend({}, R))), i.top != null && (H.top = i.top - R.top + G), i.left != null && (H.left = i.left - R.left + e), "using" in i ? i.using.call(Q, H) : V.css(H);
      }
    }, a.fn.extend({
      // offset() relates an element's border box to the document origin
      offset: function(Q) {
        if (arguments.length)
          return Q === void 0 ? this : this.each(function(e) {
            a.offset.setOffset(this, Q, e);
          });
        var i, D, w = this[0];
        if (w)
          return w.getClientRects().length ? (i = w.getBoundingClientRect(), D = w.ownerDocument.defaultView, {
            top: i.top + D.pageYOffset,
            left: i.left + D.pageXOffset
          }) : { top: 0, left: 0 };
      },
      // position() relates an element's margin box to its offset parent's padding box
      // This corresponds to the behavior of CSS absolute positioning
      position: function() {
        if (this[0]) {
          var Q, i, D, w = this[0], e = { top: 0, left: 0 };
          if (a.css(w, "position") === "fixed")
            i = w.getBoundingClientRect();
          else {
            for (i = this.offset(), D = w.ownerDocument, Q = w.offsetParent || D.documentElement; Q && (Q === D.body || Q === D.documentElement) && a.css(Q, "position") === "static"; )
              Q = Q.parentNode;
            Q && Q !== w && Q.nodeType === 1 && (e = a(Q).offset(), e.top += a.css(Q, "borderTopWidth", !0), e.left += a.css(Q, "borderLeftWidth", !0));
          }
          return {
            top: i.top - e.top - a.css(w, "marginTop", !0),
            left: i.left - e.left - a.css(w, "marginLeft", !0)
          };
        }
      },
      // This method will return documentElement in the following cases:
      // 1) For the element inside the iframe without offsetParent, this method will return
      //    documentElement of the parent window
      // 2) For the hidden or detached element
      // 3) For body or html element, i.e. in case of the html node - it will return itself
      //
      // but those exceptions were never presented as a real life use-cases
      // and might be considered as more preferable results.
      //
      // This logic, however, is not guaranteed and can change at any point in the future
      offsetParent: function() {
        return this.map(function() {
          for (var Q = this.offsetParent; Q && a.css(Q, "position") === "static"; )
            Q = Q.offsetParent;
          return Q || dg;
        });
      }
    }), a.each({ scrollLeft: "pageXOffset", scrollTop: "pageYOffset" }, function(Q, i) {
      var D = i === "pageYOffset";
      a.fn[Q] = function(w) {
        return HA(this, function(e, c, G) {
          var R;
          if (d(e) ? R = e : e.nodeType === 9 && (R = e.defaultView), G === void 0)
            return R ? R[i] : e[c];
          R ? R.scrollTo(
            D ? R.pageXOffset : G,
            D ? G : R.pageYOffset
          ) : e[c] = G;
        }, Q, w, arguments.length);
      };
    }), a.each(["top", "left"], function(Q, i) {
      a.cssHooks[i] = LB(
        K.pixelPosition,
        function(D, w) {
          if (w)
            return w = EI(D, i), WI.test(w) ? a(D).position()[i] + "px" : w;
        }
      );
    }), a.each({ Height: "height", Width: "width" }, function(Q, i) {
      a.each({
        padding: "inner" + Q,
        content: i,
        "": "outer" + Q
      }, function(D, w) {
        a.fn[w] = function(e, c) {
          var G = arguments.length && (D || typeof e != "boolean"), R = D || (e === !0 || c === !0 ? "margin" : "border");
          return HA(this, function(n, L, p) {
            var V;
            return d(n) ? w.indexOf("outer") === 0 ? n["inner" + Q] : n.document.documentElement["client" + Q] : n.nodeType === 9 ? (V = n.documentElement, Math.max(
              n.body["scroll" + Q],
              V["scroll" + Q],
              n.body["offset" + Q],
              V["offset" + Q],
              V["client" + Q]
            )) : p === void 0 ? (
              // Get width or height on the element, requesting but not forcing parseFloat
              a.css(n, L, R)
            ) : (
              // Set width or height on the element
              a.style(n, L, p, R)
            );
          }, i, G ? e : void 0, G);
        };
      });
    }), a.each([
      "ajaxStart",
      "ajaxStop",
      "ajaxComplete",
      "ajaxError",
      "ajaxSuccess",
      "ajaxSend"
    ], function(Q, i) {
      a.fn[i] = function(D) {
        return this.on(i, D);
      };
    }), a.fn.extend({
      bind: function(Q, i, D) {
        return this.on(Q, null, i, D);
      },
      unbind: function(Q, i) {
        return this.off(Q, null, i);
      },
      delegate: function(Q, i, D, w) {
        return this.on(i, Q, D, w);
      },
      undelegate: function(Q, i, D) {
        return arguments.length === 1 ? this.off(Q, "**") : this.off(i, Q || "**", D);
      },
      hover: function(Q, i) {
        return this.mouseenter(Q).mouseleave(i || Q);
      }
    }), a.each(
      "blur focus focusin focusout resize scroll click dblclick mousedown mouseup mousemove mouseover mouseout mouseenter mouseleave change select submit keydown keypress keyup contextmenu".split(" "),
      function(Q, i) {
        a.fn[i] = function(D, w) {
          return arguments.length > 0 ? this.on(i, null, D, w) : this.trigger(i);
        };
      }
    );
    var GC = /^[\s\uFEFF\xA0]+|([^\s\uFEFF\xA0])[\s\uFEFF\xA0]+$/g;
    a.proxy = function(Q, i) {
      var D, w, e;
      if (typeof i == "string" && (D = Q[i], i = Q, Q = D), !!l(Q))
        return w = E.call(arguments, 2), e = function() {
          return Q.apply(i || this, w.concat(E.call(arguments)));
        }, e.guid = Q.guid = Q.guid || a.guid++, e;
    }, a.holdReady = function(Q) {
      Q ? a.readyWait++ : a.ready(!0);
    }, a.isArray = Array.isArray, a.parseJSON = JSON.parse, a.nodeName = yA, a.isFunction = l, a.isWindow = d, a.camelCase = ag, a.type = v, a.now = Date.now, a.isNumeric = function(Q) {
      var i = a.type(Q);
      return (i === "number" || i === "string") && // parseFloat NaNs numeric-cast false positives ("")
      // ...but misinterprets leading-number strings, particularly hex literals ("0x...")
      // subtraction forces infinities to NaN
      !isNaN(Q - parseFloat(Q));
    }, a.trim = function(Q) {
      return Q == null ? "" : (Q + "").replace(GC, "$1");
    };
    var hC = A.jQuery, MC = A.$;
    return a.noConflict = function(Q) {
      return A.$ === a && (A.$ = MC), Q && A.jQuery === a && (A.jQuery = hC), a;
    }, typeof g > "u" && (A.jQuery = A.$ = a), a;
  });
})(sQ);
var qC = sQ.exports;
const O = /* @__PURE__ */ tQ(qC);
let dC = function() {
  function B(A, g) {
    this.domEl = O('<div class="aladin-popup-container"><div class="aladin-popup"><a class="aladin-closeBtn">&times;</a><div class="aladin-popupTitle"></div><div class="aladin-popupText"></div></div><div class="aladin-popup-arrow"></div></div>'), this.domEl.appendTo(A), this.view = g;
    var C = this;
    this.domEl.find(".aladin-closeBtn").click(function() {
      C.hide();
    });
  }
  return B.prototype.hide = function() {
    this.domEl.hide(), this.view.mustClearCatalog = !0, this.view.catalogForPopup.hide(), this.view.overlayForPopup.hide();
  }, B.prototype.show = function() {
    this.domEl.show();
  }, B.prototype.setTitle = function(A) {
    this.domEl.find(".aladin-popupTitle").html(A || "");
  }, B.prototype.setText = function(A) {
    this.domEl.find(".aladin-popupText").html(A || ""), this.w = this.domEl.outerWidth(), this.h = this.domEl.outerHeight();
  }, B.prototype.setSource = function(A) {
    this.source && (this.source.popup = null), A.popup = this, this.source = A, this.setPosition(A.x, A.y);
  }, B.prototype.setPosition = function(A, g) {
    var C = A - this.w / 2, I = g - this.h;
    this.source && (I += this.source.catalog.sourceSize / 2), this.domEl[0].style.left = C + "px", this.domEl[0].style.top = I + "px";
  }, B;
}(), HC = function() {
  function B() {
  }
  return B.prototype.redraw = function(A, g, C, I) {
    A.lineWidth = 1, A.strokeStyle = "rgb(150,150,220)", A.beginPath();
    for (var E, o, t = 0, s = g.length; t < s; t++)
      E = g[t], o = E.ipix, A.moveTo(E.vx[0], E.vy[0]), A.lineTo(E.vx[1], E.vy[1]), A.lineTo(E.vx[2], E.vy[2]);
    A.stroke(), A.strokeStyle = "#FFDDDD", A.beginPath();
    for (var t = 0, s = g.length; t < s; t++)
      E = g[t], o = E.ipix, A.strokeText(I + "/" + o, (E.vx[0] + E.vx[2]) / 2, (E.vy[0] + E.vy[2]) / 2);
    A.stroke();
  }, B;
}(), uA = {
  // Zenithal
  TAN: { id: 1, fov: 180 },
  /* Gnomonic projection      */
  STG: { id: 2, fov: 360 },
  /* Stereographic projection */
  SIN: { id: 3, fov: 180 },
  /* Orthographic		         */
  ZEA: { id: 4, fov: 360 },
  /* Equal-area 		         */
  FEYE: { id: 5, fov: 190 },
  AIR: { id: 6, fov: 360 },
  //AZP: {fov: 180},
  ARC: { id: 7, fov: 360 },
  NCP: { id: 8, fov: 180 },
  // Cylindrical
  MER: { id: 9, fov: 360 },
  CAR: { id: 10, fov: 360 },
  CEA: { id: 11, fov: 360 },
  CYP: { id: 12, fov: 360 },
  // Pseudo-cylindrical
  AIT: { id: 13, fov: 360 },
  PAR: { id: 14, fov: 360 },
  SFL: { id: 15, fov: 360 },
  MOL: { id: 16, fov: 360 },
  // Conic
  COD: { id: 17, fov: 360 },
  // Hybrid
  HPX: { id: 19, fov: 360 }
}, pC = [
  // Zenithal
  "SIN",
  /* Orthographic		         */
  "TAN",
  /* Gnomonic projection      */
  "STG",
  /* Stereographic projection */
  "ZEA",
  /* Equal-area 		         */
  "FEYE",
  "AIR",
  //"AZP",
  "ARC",
  "NCP",
  // Cylindrical
  "MER",
  "CAR",
  "CEA",
  "CYP",
  // Pseudo-cylindrical
  "AIT",
  "MOL",
  "PAR",
  "SFL",
  // Conic
  "COD",
  // Hybrid
  "HPX"
];
const eQ = "https://alaskybis.cds.unistra.fr/cgi/JSONProxy";
let Z = {};
Z.HTTPS_WHITELIST = [
  "alasky.u-strasbg.fr",
  "alaskybis.u-strasbg.fr",
  "alasky.unistra.fr",
  "alaskybis.unistra.fr",
  "alasky.cds.unistra.fr",
  "alaskybis.cds.unistra.fr",
  "hips.astron.nl",
  "jvo.nao.ac.jp",
  "archive.cefca.es",
  "cade.irap.omp.eu",
  "skies.esac.esa.int"
];
Z.cssScale = void 0;
Z.relMouseCoords = function(B, A) {
  if (A.offsetX)
    return { x: A.offsetX, y: A.offsetY };
  if (!Z.cssScale) {
    var g = window.getComputedStyle(document.body, null), C = g.getPropertyValue("-webkit-transform") || g.getPropertyValue("-moz-transform") || g.getPropertyValue("-ms-transform") || g.getPropertyValue("-o-transform") || g.getPropertyValue("transform"), I = /matrix\((-?\d*\.?\d+),\s*0,\s*0,\s*(-?\d*\.?\d+),\s*0,\s*0\)/, E = C.match(I);
    E ? Z.cssScale = parseFloat(E[1]) : Z.cssScale = 1;
  }
  var o = A, t = o.target || o.srcElement, s = t.currentStyle || window.getComputedStyle(t, null), M = parseInt(s.borderLeftWidth, 10), N = parseInt(s.borderTopWidth, 10), k = t.getBoundingClientRect(), U = o.clientX, Y = o.clientY;
  o.originalEvent.changedTouches && (U = o.originalEvent.changedTouches[0].clientX, Y = o.originalEvent.changedTouches[0].clientY);
  var K = U - M - k.left, l = Y - N - k.top;
  return { x: Math.round(K / Z.cssScale), y: Math.round(l / Z.cssScale) };
};
Function.prototype.bind || (Function.prototype.bind = function(B) {
  if (typeof this != "function")
    throw new TypeError("Function.prototype.bind - what is trying to be bound is not callable");
  var A = [].slice, g = A.call(arguments, 1), C = this, I = function() {
  }, E = function() {
    return C.apply(
      this instanceof I ? this : B || {},
      g.concat(A.call(arguments))
    );
  };
  return E.prototype = this.prototype, E;
});
Z.urlParam = function(B, A) {
  return A === void 0 && (A = location.search), decodeURIComponent((new RegExp("[?|&]" + B + "=([^&;]+?)(&|#|;|$)").exec(A) || [, ""])[1].replace(/\+/g, "%20")) || null;
};
Z.isNumber = function(B) {
  return !isNaN(parseFloat(B)) && isFinite(B);
};
Z.isInt = function(B) {
  return Z.isNumber(B) && Math.floor(B) === B;
};
Z.debounce = function(B, A) {
  var g = null;
  return function() {
    var C = this, I = arguments;
    clearTimeout(g), g = setTimeout(function() {
      B.apply(C, I);
    }, A);
  };
};
Z.throttle = function(B, A, g) {
  A || (A = 250);
  var C, I;
  return function() {
    var E = g || this, o = +/* @__PURE__ */ new Date(), t = arguments;
    C && o < C + A ? (clearTimeout(I), I = setTimeout(function() {
      C = o, B.apply(E, t);
    }, A)) : (C = o, B.apply(E, t));
  };
};
Z.LRUCache = function(B) {
  this._keys = [], this._items = {}, this._expires = {}, this._size = 0, this._maxsize = B || 1024;
};
Z.LRUCache.prototype = {
  set: function(B, A) {
    var g = this._keys, C = this._items, I = this._expires, E = this._size, o = this._maxsize;
    E >= o && (g.sort(function(t, s) {
      return I[t] > I[s] ? -1 : I[t] < I[s] ? 1 : 0;
    }), E--), g[E] = B, C[B] = A, I[B] = Date.now(), E++, this._keys = g, this._items = C, this._expires = I, this._size = E;
  },
  get: function(B) {
    var A = this._items[B];
    return A && (this._expires[B] = Date.now()), A;
  },
  keys: function() {
    return this._keys;
  }
};
Z.loadFromMirrors = function(B, A) {
  const g = A && A.contentType || "application/json", C = A && A.data || void 0, I = A && A.timeout || 5e3;
  if (B.length === 0)
    return Promise.reject("None of the urls given can be fetched!");
  const E = new AbortController(), o = setTimeout(() => E.abort(), I), t = {
    // *GET, POST, PUT, DELETE, etc.
    method: "GET",
    headers: {
      "Content-Type": g
    },
    // no-cors, *cors, same-origin
    mode: "cors",
    // *default, no-cache, reload, force-cache, only-if-cached
    cache: "default",
    // manual, *follow, error
    redirect: "follow",
    // Abort the request when a timeout exceeded
    signal: E.signal
  }, s = B[0] + "?" + new URLSearchParams(C);
  return fetch(s, t).then((M) => (clearTimeout(o), M.ok ? M : Promise.reject("Url: ", B[0], " cannot be reached in some way."))).catch((M) => Z.loadFromMirrors(B.slice(1), A));
};
Z.getAjaxObject = function(B, A, g, C) {
  if (C !== !1 && (C = !0), C === !0)
    var I = eQ + "?url=" + encodeURIComponent(B);
  else
    I = B;
  return A = A || "GET", g = g || null, $.ajax({
    url: I,
    method: A,
    dataType: g
  });
};
Z.isFileContext = function() {
  return window.location.protocol === "file:";
};
Z.isHttpsContext = function() {
  return window.location.protocol === "https:";
};
Z.isHttpContext = function() {
  return window.location.protocol === "http:";
};
Z.fixURLForHTTPS = function(B) {
  return Z.HTTPS_WHITELIST.some((g) => B.includes(g)) ? B.replace("http://", "https://") : B;
};
Z.getAbsoluteURL = function(B) {
  var A = document.createElement("a");
  return A.href = B, A.href;
};
Z.uuidv4 = function() {
  return "xxxxxxxx-xxxx-4xxx-yxxx-xxxxxxxxxxxx".replace(/[xy]/g, function(B) {
    var A = Math.random() * 16 | 0, g = B == "x" ? A : A & 3 | 8;
    return g.toString(16);
  });
};
Z.clone = function(B) {
  return Object.assign(
    Object.create(
      // Set the prototype of the new object to the prototype of the instance.
      // Used to allow new object behave like class instance.
      Object.getPrototypeOf(B)
    ),
    // Prevent shallow copies of nested structures like arrays, etc
    JSON.parse(JSON.stringify(B))
  );
};
Z.getDroppedFilesHandler = function(B) {
  B.preventDefault();
  let A;
  return B.dataTransfer.items ? A = [...B.dataTransfer.items] : A = [...B.dataTransfer.files], A.filter((C) => C.kind === "file").map((C) => C.getAsFile());
};
Z.dragOverHandler = function(B) {
  B.preventDefault();
};
Z.requestCORSIfNotSameOrigin = function(B) {
  return new URL(B, window.location.href).origin !== window.location.origin;
};
Z.handleCORSNotSameOrigin = function(B) {
  if (Z.requestCORSIfNotSameOrigin(B)) {
    let A = new URL(eQ);
    A.searchParams.append("url", B), B = A;
  }
  return B;
};
Z.deepCopy = function(B) {
  return Object.assign(Object.create(Object.getPrototypeOf(B)), B);
};
Z.download = function(B, A = void 0) {
  const g = document.createElement("a");
  g.href = B, A && (g.download = A), g.click();
};
let GA = function() {
};
GA.D2R = Math.PI / 180;
GA.R2D = 180 / Math.PI;
GA.sign = function(B) {
  return B > 0 ? 1 : B < 0 ? -1 : 0;
};
GA.cosd = function(B) {
  if (B % 90 == 0) {
    var A = Math.abs(Math.floor(B / 90 + 0.5)) % 4;
    switch (A) {
      case 0:
        return 1;
      case 1:
        return 0;
      case 2:
        return -1;
      case 3:
        return 0;
    }
  }
  return Math.cos(B * GA.D2R);
};
GA.sind = function(B) {
  if (B % 90 === 0) {
    var A = Math.abs(Math.floor(B / 90 - 0.5)) % 4;
    switch (A) {
      case 0:
        return 1;
      case 1:
        return 0;
      case 2:
        return -1;
      case 3:
        return 0;
    }
  }
  return Math.sin(B * GA.D2R);
};
GA.tand = function(B) {
  var A;
  return A = B % 360, A == 0 || Math.abs(A) == 180 ? 0 : A == 45 || A == 225 ? 1 : A == -135 || A == -315 ? -1 : Math.tan(B * GA.D2R);
};
GA.asind = function(B) {
  return Math.asin(B) * GA.R2D;
};
GA.acosd = function(B) {
  return Math.acos(B) * GA.R2D;
};
GA.atand = function(B) {
  return Math.atan(B) * GA.R2D;
};
GA.atan2 = function(B, A) {
  if (B != 0) {
    var g = GA.sign(B);
    if (A != 0) {
      var C = Math.atan(Math.abs(B / A));
      if (A > 0)
        return C * g;
      if (A < 0)
        return (Math.PI - C) * g;
    } else
      return Math.PI / 2 * g;
  } else
    return A > 0 ? 0 : A < 0 ? Math.PI : 0 / 0;
};
GA.atan2d = function(B, A) {
  return GA.atan2(B, A) * GA.R2D;
};
GA.cosh = function(B) {
  return (Math.exp(B) + Math.exp(-B)) / 2;
};
GA.sinh = function(B) {
  return (Math.exp(B) - Math.exp(-B)) / 2;
};
GA.tanh = function(B) {
  return (Math.exp(B) - Math.exp(-B)) / (Math.exp(B) + Math.exp(-B));
};
GA.acosh = function(B) {
  return Math.log(B + Math.sqrt(B * B - 1));
};
GA.asinh = function(B) {
  return Math.log(B + Math.sqrt(B * B + 1));
};
GA.atanh = function(B) {
  return 0.5 * Math.log((1 + B) / (1 - B));
};
GA.sinc = function(B) {
  var A = Math.abs(B), g;
  return A <= 1e-3 ? (A *= A, g = 1 - A * (1 - A / 20) / 6) : g = Math.sin(A) / A, g;
};
GA.asinc = function(B) {
  var A = Math.abs(B), g;
  return A <= 1e-3 ? (A *= A, g = 1 + A * (6 + A * (9 / 20)) / 6) : g = Math.asin(A) / A, g;
};
GA.hypot = function(B, A) {
  return Math.sqrt(B * B + A * A);
};
GA.eulerMatrix = function(B, A, g) {
  var C = new Array(3);
  C[0] = new Array(3), C[1] = new Array(3), C[2] = new Array(3);
  var I = GA.cosd(B), E = GA.sind(B), o = GA.cosd(A), t = GA.sind(A), s = GA.cosd(g), M = GA.sind(g);
  return C[0][0] = s * o * I - M * E, C[0][1] = -M * o * I - s * E, C[0][2] = -t * I, C[1][0] = s * o * E + M * I, C[1][1] = -M * o * E + s * I, C[1][2] = -t * E, C[2][0] = -t * s, C[2][1] = -t * I, C[2][2] = o, C;
};
GA.displayMatrix = function(B) {
  for (var A = B.length, g = 0, C = 0; C < A; C++)
    B[C].length > g && (g = B[C].length);
  for (var I = `<table>
`, C = 0; C < A; C++) {
    I += "<tr>";
    for (var E = 0; E < A; E++)
      I += "<td>", C < B[C].length && (I += B[C][E].toString()), I += "</td>";
    I += `</td>
`;
  }
  return I += `</table>
`, I;
};
let $A = function(B, A, g) {
  this.lon = B, this.lat = A, this.prec = g, this.frame = null, this.computeDirCos();
};
$A.factor = [3600, 60, 1];
$A.prototype = {
  setFrame: function(B) {
    this.frame = B;
  },
  computeDirCos: function() {
    var B = GA.cosd(this.lat);
    this.x = B * GA.cosd(this.lon), this.y = B * GA.sind(this.lon), this.z = GA.sind(this.lat);
  },
  computeLonLat: function() {
    var B = this.x * this.x + this.y * this.y;
    this.lon = 0, B == 0 ? this.z == 0 ? (this.lon = 0 / 0, this.lat = 0 / 0) : this.lat = this.z > 0 ? 90 : -90 : (this.lon = GA.atan2d(this.y, this.x), this.lat = GA.atan2d(this.z, Math.sqrt(B)), this.lon < 0 && (this.lon += 360));
  },
  /**
    * Squared distance between 2 points (= 4.sin<sup>2</sup>(r/2))
    * @param  pos      another position on the sphere
    * @return ||pos-this||<sup>2</sup> = 4.sin<sup>2</sup>(r/2)
   **/
  dist2: function(B) {
    var A = B.x - this.x, g = A * A;
    return A = B.y - this.y, g += A * A, A = B.z - this.z, g += A * A, g;
  },
  /**
   * Distance between 2 points on the sphere.
   * @param  pos another position on the sphere
   * @return distance in degrees in range [0, 180]
  **/
  distance: function(B) {
    return B.x == 0 && B.y == 0 && B.z == 0 || this.x == 0 && this.y == 0 && this.z == 0 ? 0 / 0 : 2 * GA.asind(0.5 * Math.sqrt(this.dist2(B)));
  },
  /**
   * Transform the position into another frame.
   * @param new_frame	The frame of the resulting position.
  **/
  convertTo: function(B) {
    this.frame.equals(B) || (this.frame.toICRS(this.coo), B.fromICRS(this.coo), this.frame = B, this.lon = this.lat = 0 / 0);
  },
  /**
   * Rotate a coordinate (apply a rotation to the position).
   * @param R [3][3] Rotation Matrix
   */
  rotate: function(B) {
    var A, g, C;
    B != Umatrix3 && (A = B[0][0] * this.x + B[0][1] * this.y + B[0][2] * this.z, g = B[1][0] * this.x + B[1][1] * this.y + B[1][2] * this.z, C = B[2][0] * this.x + B[2][1] * this.y + B[2][2] * this.z, this.x = A, this.y = g, this.z = C, this.lon = this.lat = 0 / 0);
  },
  /**
   * Rotate a coordinate (apply a rotation to the position) in reverse direction.
   * The method is the inverse of rotate.
   * @param R [3][3] Rotation Matrix
   */
  rotate_1: function(B) {
    var A, g, C;
    B != Umatrix3 && (A = B[0][0] * this.x + B[1][0] * this.y + B[2][0] * this.z, g = B[0][1] * this.x + B[1][1] * this.y + B[2][1] * this.z, C = B[0][2] * this.x + B[1][2] * this.y + B[2][2] * this.z, this.x = A, this.y = g, this.z = C, this.lon = this.lat = 0 / 0);
  },
  /**
   * Test equality of Coo.
   * @param coo Second coordinate to compare with
   * @return  True if the two coordinates are equal
   */
  equals: function(B) {
    return this.x == B.x && this.y == B.y && this.z == B.z;
  },
  /**
   * parse a coordinate string. The coordinates can be in decimal or sexagesimal
   * @param str string to parse
   * @return true if the parsing succeded, false otherwise
   */
  parse: function(B) {
    var A = B.indexOf("+");
    if (A < 0 && (A = B.indexOf("-")), A < 0 && (A = B.indexOf(" ")), A < 0)
      return this.lon = 0 / 0, this.lat = 0 / 0, this.prec = 0, !1;
    var g = B.substring(0, A), C = B.substring(A);
    return this.lon = this.parseLon(g), this.lat = this.parseLat(C), !0;
  },
  parseLon: function(A) {
    var A = A.trim();
    if (A = A.replace(/:/g, " "), A.indexOf(" ") < 0) {
      var g = A.indexOf(".");
      return this.prec = g < 0 ? 0 : A.length - g - 1, parseFloat(A);
    } else {
      for (var C = new tB(A, " "), I = 0, E = 0, o = 0; C.hasMore(); ) {
        var t = C.nextToken(), s = t.indexOf(".");
        switch (E += parseFloat(t) * $A.factor[I], I) {
          case 0:
            o = s < 0 ? 1 : 2;
            break;
          case 1:
            o = s < 0 ? 3 : 4;
            break;
          case 2:
            o = s < 0 ? 5 : 4 + t.length - s;
        }
        I++;
      }
      return this.prec = o, E * 15 / 3600;
    }
  },
  parseLat: function(A) {
    var A = A.trim();
    A = A.replace(/:/g, " ");
    var g;
    if (A.charAt(0) == "-" ? (g = -1, A = A.substring(1)) : A.charAt(0) == "-" ? (g = 1, A = A.substring(1)) : g = 1, A.indexOf(" ") < 0) {
      var C = A.indexOf(".");
      return this.prec = C < 0 ? 0 : A.length - C - 1, parseFloat(A) * g;
    } else {
      for (var I = new tB(A, " "), E = 0, o = 0, t = 0; I.hasMore(); ) {
        var s = I.nextToken(), M = s.indexOf(".");
        switch (o += parseFloat(s) * $A.factor[E], E) {
          case 0:
            t = M < 0 ? 1 : 2;
            break;
          case 1:
            t = M < 0 ? 3 : 4;
            break;
          case 2:
            t = M < 0 ? 5 : 4 + s.length - M;
        }
        E++;
      }
      return this.prec = t, o * g / 3600;
    }
  },
  /**
   * Format coordinates according to the options
   * @param options 'd': decimal, 's': sexagésimal, '/': space separated, '2': return [ra,dec] in an array
   * @return the formatted coordinates
   */
  format: function(B) {
    isNaN(this.lon) && this.computeLonLat();
    var A = "", g = "";
    if (B.indexOf("d") >= 0)
      A = cg.format(this.lon, this.prec), g = cg.format(this.lat, this.prec), A.indexOf("0") == 0 && (A = A.substr(1)), g.indexOf("0") == 0 && (g = g.substr(1));
    else
      var C = this.lon / 15, A = cg.toSexagesimal(C, this.prec + 1, !1), g = cg.toSexagesimal(this.lat, this.prec, !1);
    return this.lat >= 0 && (g = "+" + g), B.indexOf("/") >= 0 ? A + " " + g : B.indexOf("2") >= 0 ? [A, g] : A + g;
  }
};
function tB(B, A) {
  this.string = cQ.trim(B, A), this.sep = A, this.pos = 0;
}
tB.prototype = {
  /**
   * Check if the string has more tokens
   * @return true if a token remains (read with nextToken())
   */
  hasMore: function() {
    return this.pos < this.string.length;
  },
  /**
   * Returns the next token (as long as hasMore() is true)
   * @return the token string
   */
  nextToken: function() {
    for (var B = this.pos; B < this.string.length && this.string.charAt(B) == this.sep; )
      B++;
    for (var A = B; A < this.string.length && this.string.charAt(A) != this.sep; )
      A++;
    return this.pos = A, this.string.substring(B, A);
  }
};
function cQ() {
}
cQ.trim = function(B, A) {
  for (var g = 0, C = B.length - 1; g < B.length && B.charAt(g) == A; )
    g++;
  if (g == B.length)
    return "";
  for (; C > g && B.charAt(C) == A; )
    C--;
  return B.substring(g, C + 1);
};
function cg() {
}
cg.pow10 = [
  1,
  10,
  100,
  1e3,
  1e4,
  1e5,
  1e6,
  1e7,
  1e8,
  1e9,
  //      10           11            12             13              14
  1e10,
  1e11,
  1e12,
  1e13,
  1e14
];
cg.rndval = [
  0.5,
  0.05,
  5e-3,
  5e-4,
  5e-5,
  5e-6,
  5e-7,
  5e-8,
  //      8            9             10             11              12
  5e-9,
  5e-10,
  5e-11,
  5e-12,
  5e-13,
  //      13                14
  5e-14,
  5e-14
];
cg.format = function(B, A) {
  if (A <= 0)
    return Math.round(B).toString();
  var g = B.toString(), C = g.indexOf("."), I = C >= 0 ? g.length - C - 1 : 0;
  if (A >= I) {
    C < 0 && (g += "."), g.length == 2 && (g = "0" + g);
    for (var E = 0; E < A - I; E++)
      g += "0";
    return g;
  }
  return g = (B + cg.rndval[A]).toString(), C == 1 && (g = "0" + g, C = C + 1), g.substr(0, C + A + 1);
};
cg.toSexagesimal = function(B, A, g) {
  var C = B < 0 ? "-" : g ? "+" : "", I = Math.abs(B);
  switch (A) {
    case 1:
      var E = Math.round(I);
      return C + E.toString();
    case 2:
      return C + cg.format(I, 1);
    case 3:
      var E = Math.floor(I), o = Math.round((I - E) * 60);
      return C + E + " " + o;
    case 4:
      var E = Math.floor(I), o = (I - E) * 60;
      return C + E + " " + cg.format(o, 1);
    case 5:
      var E = Math.floor(I), o = (I - E) * 60, t = Math.floor(o), s = Math.round((o - t) * 60);
      return C + E + " " + t + " " + s;
    case 6:
    case 7:
    case 8:
      var E = Math.floor(I);
      E < 10 && (E = "0" + E);
      var o = (I - E) * 60, t = Math.floor(o);
      t < 10 && (t = "0" + t);
      var s = (o - t) * 60;
      return C + E + " " + t + " " + cg.format(s, A - 5);
    default:
      return C + cg.format(I, 1);
  }
};
let RA = function() {
  var B = { J2000: "J2000", GAL: "Galactic" };
  return {
    SYSTEMS: B,
    J2000: { label: "J2000", system: B.J2000 },
    J2000d: { label: "J2000d", system: B.J2000 },
    GAL: { label: "Galactic", system: B.GAL },
    fromString: function(A, g) {
      return A ? (A = A.toLowerCase().replace(/^\s+|\s+$/g, ""), A.indexOf("j2000d") == 0 || A.indexOf("icrsd") == 0 ? RA.J2000d : A.indexOf("j2000") == 0 || A.indexOf("icrs") == 0 ? RA.J2000 : A.indexOf("gal") == 0 ? RA.GAL : g || null) : g || null;
    }
  };
}(), jA = function() {
  return {
    /**
     * passage de xy projection à xy dans la vue écran 
     * @param x
     * @param y
     * @param width
     * @param height
     * @param largestDim largest dimension of the view
     * @returns position in the view
     */
    xyToView: function(B, A, g, C, I, E, o) {
      return o == null && (o = !0), o ? { vx: jA.myRound(I / 2 * (1 + E * B) - (I - g) / 2), vy: jA.myRound(I / 2 * (1 + E * A) - (I - C) / 2) } : { vx: I / 2 * (1 + E * B) - (I - g) / 2, vy: I / 2 * (1 + E * A) - (I - C) / 2 };
    },
    /**
     * passage de xy dans la vue écran à xy projection
     * @param vx
     * @param vy
     * @param width
     * @param height
     * @param largestDim
     * @param zoomFactor
     * @returns position in xy projection
     */
    viewToXy: function(B, A, g, C, I, E) {
      return { x: ((2 * B + (I - g)) / I - 1) / E, y: ((2 * A + (I - C)) / I - 1) / E };
    },
    /**
     * convert a 
     * @returns position x,y in the view. Null if projection is impossible
     */
    /*radecToViewXy: function(ra, dec, currentProjection, currentFrame, width, height, largestDim, zoomFactor) {
                var xy;
                if (currentFrame.system != CooFrameEnum.SYSTEMS.J2000) {
                    var lonlat = CooConversion.J2000ToGalactic([ra, dec]);
                    xy = currentProjection.project(lonlat[0], lonlat[1]);
                }
                else {
                    xy = currentProjection.project(ra, dec);
                }
                if (!xy) {
                    return null;
                }
    
                return AladinUtils.xyToView(xy.X, xy.Y, width, height, largestDim, zoomFactor, false);
            },*/
    radecToViewXy: function(B, A, g) {
      return g.wasm.worldToScreen(B, A);
    },
    viewXyToClipXy: function(B, A, g) {
      return g.wasm.screenToClip(B, A);
    },
    myRound: function(B) {
      return B < 0 ? -1 * (-B | 0) : B | 0;
    },
    /**
     * Test whether a xy position is the view
     * @param vx
     * @param vy
     * @param width
     * @param height
     * @returns a boolean whether (vx, vy) is in the screen
     */
    isInsideViewXy: function(B, A, g, C) {
      return B >= 0 && B < g && A >= 0 && A < C;
    },
    /**
     * tests whether a healpix pixel is visible or not
     * @param pixCorners array of position (xy view) of the corners of the pixel
     * @param viewW
     */
    isHpxPixVisible: function(B, A, g) {
      for (var C = 0; C < B.length; C++)
        if (B[C].vx >= -20 && B[C].vx < A + 20 && B[C].vy >= -20 && B[C].vy < g + 20)
          return !0;
      return !1;
    },
    ipixToIpix: function(B, A, g) {
    },
    // Zoom is handled in the backend
    /*getZoomFactorForAngle: function(angleInDegrees, projectionMethod) {
                var p1 = {ra: 0, dec: 0};
                var p2 = {ra: angleInDegrees, dec: 0};
                var projection = new Projection(angleInDegrees/2, 0);
                projection.setProjection(projectionMethod);
                var p1Projected = projection.project(p1.ra, p1.dec);
                var p2Projected = projection.project(p2.ra, p2.dec);
               
                var zoomFactor = 1/Math.abs(p1Projected.X - p2Projected.Y);
    
                return zoomFactor;
            },*/
    counterClockwiseTriangle: function(B, A, g, C, I, E) {
      return B * C + A * I + g * E - I * C - E * B - g * A >= 0;
    },
    // grow array b of vx,vy view positions by *val* pixels
    grow2: function(B, A) {
      for (var g = 0, C = 0; C < 4; C++)
        B[C] == null && g++;
      if (g > 1)
        return B;
      for (var I = [], C = 0; C < 4; C++)
        I.push({ vx: B[C].vx, vy: B[C].vy });
      for (var C = 0; C < 2; C++) {
        var E = C == 1 ? 1 : 0, o = C == 1 ? 3 : 2;
        if (I[E] == null) {
          var t, s;
          E == 0 || E == 3 ? (t = 1, s = 2) : (t = 0, s = 3), I[E] = { vx: (I[t].vx + I[s].vx) / 2, vy: (I[t].vy + I[s].vy) / 2 };
        }
        if (I[o] == null) {
          var t, s;
          o == 0 || o == 3 ? (t = 1, s = 2) : (t = 0, s = 3), I[o] = { vx: (I[t].vx + I[s].vx) / 2, vy: (I[t].vy + I[s].vy) / 2 };
        }
        if (!(I[E] == null || I[o] == null)) {
          var M = Math.atan2(I[o].vy - I[E].vy, I[o].vx - I[E].vx), N = A * Math.cos(M);
          I[E].vx -= N, I[o].vx += N, N = A * Math.sin(M), I[E].vy -= N, I[o].vy += N;
        }
      }
      return I;
    },
    // SVG icons templates are stored here rather than in a CSS, as to allow
    // to dynamically change the fill color
    // Pretty ugly, haven't found a prettier solution yet
    //
    // TODO: store this in the Stack class once it will exist
    //
    SVG_ICONS: {
      CATALOG: '<svg xmlns="http://www.w3.org/2000/svg"><polygon points="1,0,5,0,5,3,1,3"  fill="FILLCOLOR" /><polygon points="7,0,9,0,9,3,7,3"  fill="FILLCOLOR" /><polygon points="10,0,12,0,12,3,10,3"  fill="FILLCOLOR" /><polygon points="13,0,15,0,15,3,13,3"  fill="FILLCOLOR" /><polyline points="1,5,5,9"  stroke="FILLCOLOR" /><polyline points="1,9,5,5" stroke="FILLCOLOR" /><line x1="7" y1="7" x2="15" y2="7" stroke="FILLCOLOR" stroke-width="2" /><polyline points="1,11,5,15"  stroke="FILLCOLOR" /><polyline points="1,15,5,11"  stroke="FILLCOLOR" /><line x1="7" y1="13" x2="15" y2="13" stroke="FILLCOLOR" stroke-width="2" /></svg>',
      MOC: '<svg xmlns="http://www.w3.org/2000/svg"><polyline points="0.5,7,2.5,7,2.5,5,7,5,7,3,10,3,10,5,13,5,13,7,15,7,15,9,13,9,13,12,10,12,10,14,7,14,7,12,2.5,12,2.5,10,0.5,10,0.5,7" stroke-width="1" stroke="FILLCOLOR" fill="transparent" /><line x1="1" y1="10" x2="6" y2="5" stroke="FILLCOLOR" stroke-width="0.5" /><line x1="2" y1="12" x2="10" y2="4" stroke="FILLCOLOR" stroke-width="0.5" /><line x1="5" y1="12" x2="12" y2="5" stroke="FILLCOLOR" stroke-width="0.5" /><line x1="7" y1="13" x2="13" y2="7" stroke="FILLCOLOR" stroke-width="0.5" /><line x1="10" y1="13" x2="13" y2="10" stroke="FILLCOLOR" stroke-width="0.5" /></svg>',
      OVERLAY: '<svg xmlns="http://www.w3.org/2000/svg"><polygon points="10,5,10,1,14,1,14,14,2,14,2,9,6,9,6,5" fill="transparent" stroke="FILLCOLOR" stroke-width="2"/></svg>'
    },
    /**
    * @function degreesToString
    * Convert a number in degrees into a string<br>
    *
    * @param numberDegrees number in degrees (integer or decimal)
    * @return a formattes string
    * 
    * @example <caption> Result in degrees </caption>
    * // returns "1°"
    * Numbers.degreesToString(1)
    * @example <caption> Result in arcminutes </caption>
    * // returns "6 arcmin"
    * Numbers.degreesToString(0.1);
    * @example <caption> Result in arcseconds </caption>
    * // returns "36 arcsec"
    * Numbers.degreesToString(0.01);
    */
    degreesToString: function(B) {
      let A = 3, g = B | 0, C = Math.abs(B - g), I = C * 60 | 0, E = ((C * 60 - I) * 60).toPrecision(A);
      return g != 0 ? B.toPrecision(A) + "°" : I != 0 ? (C * 60).toPrecision(A) + " arcmin" : E != 0 ? E + " arcsec" : "0°";
    }
  };
}(), fC = function() {
  const B = {};
  return B.MIRRORS = ["https://alasky.cds.unistra.fr/cgi/simbad-flat/simbad-quick.py", "https://alaskybis.cds.unistra.fr/cgi/simbad-flat/simbad-quick.py"], B.query = function(A, g, C, I) {
    var E = new $A(A, g, 7), o = { Ident: E.format("s/"), SR: C };
    Z.loadFromMirrors(B.MIRRORS, { contentType: "text/plain", data: o }).then((t) => t.text()).then((t) => {
      I.view.setCursor("pointer");
      var s = /(.*?)\/(.*?)\((.*?),(.*?)\)/g, M = s.exec(t);
      if (M) {
        var N = new $A();
        N.parse(M[1]);
        var k = M[2], U = '<div class="aladin-sp-title"><a target="_blank" href="https://simbad.cds.unistra.fr/simbad/sim-id?Ident=' + encodeURIComponent(k) + '">' + k + "</a></div>", Y = '<div class="aladin-sp-content">';
        Y += "<em>Type: </em>" + M[4] + "<br>";
        var K = M[3];
        Z.isNumber(K) && (Y += "<em>Mag: </em>" + K + "<br>"), Y += '<br><a target="_blank" href="http://cdsportal.u-strasbg.fr/?target=' + encodeURIComponent(k) + '">Query in CDS portal</a>', Y += "</div>", I.showPopup(N.lon, N.lat, U, Y);
      } else {
        let l = '<div class="aladin-sp-title">Ohoh</div>', q = '<div class="aladin-sp-content">No match was found on <a href="https://simbad.cds.unistra.fr/simbad">Simbad</a> in ' + jA.degreesToString(C) + " around this point.</div>";
        I.showPopup(E.lon, E.lat, l, q);
      }
    }).catch(
      (t) => {
        I.view.setCursor("pointer"), I.hidePopup();
      }
    );
  }, B;
}(), uC = function() {
  const B = {};
  return B.MIRRORS = ["https://alasky.cds.unistra.fr/planetary-features/cs"], B.PLANETS_RADIUS = {
    mercury: 2439400,
    venus: 6051e3,
    earth: 6378137,
    mars: 3396190,
    moon: 1738100,
    ceres: 473e3,
    titan: 2575e3,
    titania: 788400,
    dione: 561400,
    enceladus: 252100,
    iapetus: 734500,
    mimas: 198200,
    rhea: 763800,
    tethys: 533e3,
    callisto: 2410300,
    ariel: 578900,
    charon: 606e3,
    triton: 1353e3,
    pluto: 1188300
  }, B.HAS_WEST_INCREASING_LONGITUDES = ["amalthea", "callisto", "deimos", "dione", "enceladus", "epimetheus", "eros", "europa", "ganymede", "gaspra", "hyperion", "iapetus", "io", "janus", "mathilde", "mercury", "mimas", "phobos", "phoebe", "proteus", "rhea", "tethys", "thebe", "titan"], B.query = function(A, g, C, I, E) {
    let o = A;
    const t = g;
    B.HAS_WEST_INCREASING_LONGITUDES.includes(I) && (o = 360 - o);
    const s = { lon: o, lat: t, SR: C, format: "csv", body: I };
    function M(N) {
      let k = "", U = [""], Y = [U], K = 0, l = 0, d = !0, q;
      for (q of N)
        q === '"' ? (d && q === k && (U[K] += q), d = !d) : q === "," && d ? q = U[++K] = "" : q === `
` && d ? (k === "\r" && (U[K] = U[K].slice(0, -1)), U = Y[++l] = [q = ""], K = 0) : U[K] += q, k = q;
      return Y;
    }
    Z.loadFromMirrors(B.MIRRORS, { contentType: "text/plain", data: s }).then((N) => N.text()).then((N) => {
      E.view.setCursor("pointer");
      const k = N.split(`
`), U = M(k[0])[0];
      if (k.length > 1 && k[1].length > 0) {
        const Y = M(k[1])[0], K = U.findIndex((hA) => hA.includes("longitude")), l = U.findIndex((hA) => hA.includes("latitude")), d = U.findIndex((hA) => hA.includes("feature_name")), q = Y[d], W = Y[U.indexOf("feature_id")], T = '<div class="aladin-sp-title"><a target="_blank" href="https://planetarynames.wr.usgs.gov/Feature/' + W + '">' + q + "</a></div>", v = Y[U.indexOf("feature_type")];
        let _ = '<div class="aladin-sp-content"> </div>';
        _ += "<em>Type: </em>" + v + "<br><br>", _ += '<a target="_blank" href="https://planetarynames.wr.usgs.gov/Feature/' + W + '">More information</a>';
        let a = parseFloat(Y[K]);
        const iA = parseFloat(Y[l]);
        let IA = !0;
        const z = U.indexOf("coordinate_system");
        z > 0 && Y[z].includes("+West") && (IA = !1), IA || (a = 360 - a);
        let oA;
        try {
          const hA = parseFloat(Y[U.indexOf("diameter")]);
          if (I in B.PLANETS_RADIUS) {
            const yA = 2 * Math.PI * B.PLANETS_RADIUS[I] * Math.cos(iA * Math.PI / 8180);
            oA = 180 * (2 * Math.PI * (1e3 * hA / 2) / yA) / Math.PI;
          }
        } catch (hA) {
          console.error(hA);
        }
        E.showPopup(a, iA, T, _, oA);
      } else {
        let Y = '<div class="aladin-sp-title">Ohoh</div>', l = '<div class="aladin-sp-content">No match was found on <a href="https://planetarynames.wr.usgs.gov">planetarynames.wr.usgs.gov</a> in ' + jA.degreesToString(C) + " around this point.</div>";
        E.showPopup(A, g, Y, l);
      }
    }).catch((N) => {
      E.view.setCursor("pointer"), E.hidePopup();
    });
  }, B;
}(), GQ = function(B, A) {
  const g = Z.relMouseCoords(B.imageCanvas, A);
  let C = B.wasm.screenToWorld(g.x, g.y);
  if (C)
    if (B.aladin.getBaseImageLayer().isPlanetaryBody() === !1) {
      const I = Math.min(1, 15 * B.fov / B.largestDim);
      console.log('queryRadius "generic pointer": ', I), fC.query(C[0], C[1], I, B.aladin);
    } else {
      const I = B.aladin.getBaseImageLayer().properties.hipsBody;
      uC.query(C[0], C[1], Math.min(80, B.fov / 20), I, B.aladin);
    }
  else
    console.log("The location you clicked on is out of the view.");
}, xC = function() {
  function B(nA, og, PA) {
    var gg, XA, OA;
    for (XA = 0; XA < 30; XA++)
      for (gg = 0; gg < 73; gg++)
        OA = (gg + XA * 74) * 4, nA[OA] = nA[OA + 4], nA[OA + 1] = nA[OA + 5], nA[OA + 2] = nA[OA + 6];
    for (XA = 0; XA < 30; XA++)
      OA = (73 + XA * 74) * 4, XA < og ? (nA[OA] = gA[PA].bg.r, nA[OA + 1] = gA[PA].bg.g, nA[OA + 2] = gA[PA].bg.b) : (nA[OA] = gA[PA].fg.r, nA[OA + 1] = gA[PA].fg.g, nA[OA + 2] = gA[PA].fg.b);
  }
  var A = 0, g = 2, C, I = 0, E = (/* @__PURE__ */ new Date()).getTime(), o = E, t = E, s = 0, M = 1e3, N = 0, k, U, Y, K, l, d = 0, q = 1e3, W = 0, T, v, _, a, iA = 0, IA = 1e3, z = 0, oA, hA, yA, rA, gA = { fps: { bg: { r: 16, g: 16, b: 48 }, fg: { r: 0, g: 255, b: 255 } }, ms: { bg: { r: 16, g: 48, b: 16 }, fg: { r: 0, g: 255, b: 0 } }, mb: { bg: {
    r: 48,
    g: 16,
    b: 26
  }, fg: { r: 255, g: 0, b: 128 } } };
  C = document.createElement("div"), C.style.cursor = "pointer", C.style.width = "80px", C.style.opacity = "0.9", C.style.zIndex = "10001", C.addEventListener("click", function() {
    switch (A++, A == g && (A = 0), k.style.display = "none", T.style.display = "none", oA.style.display = "none", A) {
      case 0:
        k.style.display = "block";
        break;
      case 1:
        T.style.display = "block";
        break;
      case 2:
        oA.style.display = "block";
    }
  }, !1), k = document.createElement("div"), k.style.backgroundColor = "rgb(" + Math.floor(gA.fps.bg.r / 2) + "," + Math.floor(gA.fps.bg.g / 2) + "," + Math.floor(gA.fps.bg.b / 2) + ")", k.style.padding = "2px 0px 3px 0px", C.appendChild(k), U = document.createElement("div"), U.style.fontFamily = "Helvetica, Arial, sans-serif", U.style.textAlign = "left", U.style.fontSize = "9px", U.style.color = "rgb(" + gA.fps.fg.r + "," + gA.fps.fg.g + "," + gA.fps.fg.b + ")", U.style.margin = "0px 0px 1px 3px", U.innerHTML = '<span style="font-weight:bold">FPS</span>', k.appendChild(U), Y = document.createElement("canvas"), Y.width = 74, Y.height = 30, Y.style.display = "block", Y.style.marginLeft = "3px", k.appendChild(Y), K = Y.getContext("2d"), K.fillStyle = "rgb(" + gA.fps.bg.r + "," + gA.fps.bg.g + "," + gA.fps.bg.b + ")", K.fillRect(0, 0, Y.width, Y.height), l = K.getImageData(0, 0, Y.width, Y.height), T = document.createElement("div"), T.style.backgroundColor = "rgb(" + Math.floor(gA.ms.bg.r / 2) + "," + Math.floor(gA.ms.bg.g / 2) + "," + Math.floor(gA.ms.bg.b / 2) + ")", T.style.padding = "2px 0px 3px 0px", T.style.display = "none", C.appendChild(T), v = document.createElement("div"), v.style.fontFamily = "Helvetica, Arial, sans-serif", v.style.textAlign = "left", v.style.fontSize = "9px", v.style.color = "rgb(" + gA.ms.fg.r + "," + gA.ms.fg.g + "," + gA.ms.fg.b + ")", v.style.margin = "0px 0px 1px 3px", v.innerHTML = '<span style="font-weight:bold">MS</span>', T.appendChild(v), Y = document.createElement("canvas"), Y.width = 74, Y.height = 30, Y.style.display = "block", Y.style.marginLeft = "3px", T.appendChild(Y), _ = Y.getContext("2d"), _.fillStyle = "rgb(" + gA.ms.bg.r + "," + gA.ms.bg.g + "," + gA.ms.bg.b + ")", _.fillRect(0, 0, Y.width, Y.height), a = _.getImageData(0, 0, Y.width, Y.height);
  try {
    performance && performance.memory && performance.memory.totalJSHeapSize && (g = 3);
  } catch {
  }
  return oA = document.createElement("div"), oA.style.backgroundColor = "rgb(" + Math.floor(gA.mb.bg.r / 2) + "," + Math.floor(gA.mb.bg.g / 2) + "," + Math.floor(gA.mb.bg.b / 2) + ")", oA.style.padding = "2px 0px 3px 0px", oA.style.display = "none", C.appendChild(oA), hA = document.createElement("div"), hA.style.fontFamily = "Helvetica, Arial, sans-serif", hA.style.textAlign = "left", hA.style.fontSize = "9px", hA.style.color = "rgb(" + gA.mb.fg.r + "," + gA.mb.fg.g + "," + gA.mb.fg.b + ")", hA.style.margin = "0px 0px 1px 3px", hA.innerHTML = '<span style="font-weight:bold">MB</span>', oA.appendChild(hA), Y = document.createElement("canvas"), Y.width = 74, Y.height = 30, Y.style.display = "block", Y.style.marginLeft = "3px", oA.appendChild(Y), yA = Y.getContext("2d"), yA.fillStyle = "#301010", yA.fillRect(0, 0, Y.width, Y.height), rA = yA.getImageData(0, 0, Y.width, Y.height), { domElement: C, update: function() {
    I++, E = (/* @__PURE__ */ new Date()).getTime(), d = E - o, q = Math.min(q, d), W = Math.max(W, d), B(a.data, Math.min(30, 30 - d / 200 * 30), "ms"), v.innerHTML = '<span style="font-weight:bold">' + d + " MS</span> (" + q + "-" + W + ")", _.putImageData(a, 0, 0), o = E, E > t + 1e3 && (s = Math.round(I * 1e3 / (E - t)), M = Math.min(M, s), N = Math.max(N, s), B(l.data, Math.min(30, 30 - s / 100 * 30), "fps"), U.innerHTML = '<span style="font-weight:bold">' + s + " FPS</span> (" + M + "-" + N + ")", K.putImageData(l, 0, 0), g == 3 && (iA = performance.memory.usedJSHeapSize * 954e-9, IA = Math.min(IA, iA), z = Math.max(z, iA), B(rA.data, Math.min(30, 30 - iA / 2), "mb"), hA.innerHTML = '<span style="font-weight:bold">' + Math.round(iA) + " MB</span> (" + Math.round(IA) + "-" + Math.round(z) + ")", yA.putImageData(rA, 0, 0)), t = E, I = 0);
  } };
}, gI = function() {
  let B = function(A) {
    A = A || {}, this.uuid = Z.uuidv4(), this.type = "overlay", this.name = A.name || "overlay", this.color = A.color || Color.getNextColor(), this.lineWidth = A.lineWidth || 2, this.overlayItems = [], this.isShowing = !0;
  };
  return B.prototype.show = function() {
    this.isShowing || (this.isShowing = !0, this.overlayItems.forEach((A) => A.show()), this.reportChange());
  }, B.prototype.hide = function() {
    this.isShowing && (this.isShowing = !1, this.overlayItems.forEach((A) => A.hide()), this.reportChange());
  }, B.parseSTCS = function(A, g) {
    g = g || {};
    for (var C = [], I = A.match(/\S+/g), E = 0, o = I.length; E < o; ) {
      var t = I[E].toLowerCase();
      if (t == "polygon") {
        var s = [];
        if (E++, k = I[E].toLowerCase(), Z.isNumber(k) && (k = "icrs", E--), k == "icrs" || k == "j2000" || k == "fk5") {
          for (; E + 2 < o; ) {
            var M = parseFloat(I[E + 1]);
            if (isNaN(M))
              break;
            var N = parseFloat(I[E + 2]);
            s.push([M, N]), E += 2;
          }
          g.closed = !0, C.push(zA.polygon(s, g));
        }
      } else if (t == "circle") {
        var k;
        if (E++, k = I[E].toLowerCase(), Z.isNumber(k) && (k = "icrs", E--), k == "icrs" || k == "j2000" || k == "fk5") {
          var M, N, U;
          M = parseFloat(I[E + 1]), N = parseFloat(I[E + 2]), U = parseFloat(I[E + 3]), C.push(zA.circle(M, N, U, g)), E += 3;
        }
      }
      E++;
    }
    return C;
  }, B.prototype.addFootprints = function(A) {
    for (var g = 0, C = A.length; g < C; g++)
      this.add(A[g], !1);
  }, B.prototype.add = function(A, g) {
    g = g !== void 0 ? g : !0, this.overlayItems.push(A), A.setOverlay(this), g && this.view.requestRedraw();
  }, B.prototype.getFootprint = function(A) {
    return A < this.footprints.length ? this.footprints[A] : null;
  }, B.prototype.setView = function(A) {
    this.view = A;
  }, B.prototype.removeAll = function() {
    this.overlayItems = [];
  }, B.prototype.draw = function(A) {
    if (this.isShowing) {
      A.strokeStyle = this.color, A.lineWidth = this.lineWidth;
      for (var g = 0; g < this.overlayItems.length; g++)
        this.overlayItems[g].draw(A, this.view);
    }
  }, B.increaseBrightness = function(A, g) {
    A = A.replace(/^\s*#|\s*$/g, ""), A.length == 3 && (A = A.replace(/(.)/g, "$1$1"));
    var C = parseInt(A.substr(0, 2), 16), I = parseInt(A.substr(2, 2), 16), E = parseInt(A.substr(4, 2), 16);
    return "#" + (0 | 256 + C + (256 - C) * g / 100).toString(16).substr(1) + (0 | 256 + I + (256 - I) * g / 100).toString(16).substr(1) + (0 | 256 + E + (256 - E) * g / 100).toString(16).substr(1);
  }, B.prototype.setColor = function(A) {
    this.color = A, this.reportChange();
  }, B.prototype.setLineWidth = function(A) {
    this.lineWidth = A, this.reportChange();
  }, B.prototype.reportChange = function() {
    this.view && this.view.requestRedraw();
  }, B;
}(), hQ = function() {
  let B = function(A, g, C) {
    C = C || {}, this.color = C.color || void 0, this.fillColor = C.fillColor || void 0, this.lineWidth = C.lineWidth || 2, this.id = "circle-" + Z.uuidv4(), this.setCenter(A), this.setRadius(g), this.overlay = null, this.isShowing = !0, this.isSelected = !1, this.selectionColor = "#00ff00";
  };
  return B.prototype.setColor = function(A) {
    this.color != A && (this.color = A, this.overlay && this.overlay.reportChange());
  }, B.prototype.setSelectionColor = function(A) {
    this.selectionColor != A && (this.selectionColor = A, this.overlay && this.overlay.reportChange());
  }, B.prototype.setLineWidth = function(A) {
    this.lineWidth != A && (this.lineWidth = A, this.overlay && this.overlay.reportChange());
  }, B.prototype.setOverlay = function(A) {
    this.overlay = A;
  }, B.prototype.show = function() {
    this.isShowing || (this.isShowing = !0, this.overlay && this.overlay.reportChange());
  }, B.prototype.hide = function() {
    this.isShowing && (this.isShowing = !1, this.overlay && this.overlay.reportChange());
  }, B.prototype.select = function() {
    this.isSelected || (this.isSelected = !0, this.overlay && this.overlay.reportChange());
  }, B.prototype.deselect = function() {
    this.isSelected && (this.isSelected = !1, this.overlay && this.overlay.reportChange());
  }, B.prototype.isFootprint = function() {
    return !0;
  }, B.prototype.setCenter = function(A) {
    this.centerRaDec = A, this.overlay && this.overlay.reportChange();
  }, B.prototype.setRadius = function(A) {
    this.radiusDegrees = A, this.overlay && this.overlay.reportChange();
  }, B.prototype.draw = function(A, g, C) {
    if (!this.isShowing)
      return;
    C = C === !0 || !1;
    var I = jA.radecToViewXy(this.centerRaDec[0], this.centerRaDec[1], g);
    if (!I)
      return;
    this.center = {
      x: I[0],
      y: I[1]
    };
    var E = this.centerRaDec[0], o = this.centerRaDec[1] + (E > 0 ? -this.radiusDegrees : this.radiusDegrees);
    let t = jA.radecToViewXy(E, o, g);
    if (!t)
      return;
    let [s, M] = jA.viewXyToClipXy(this.center.x, this.center.y, g), [N, k] = jA.viewXyToClipXy(t[0], t[1], g);
    if (!((s - N) * (s - N) + (M - k) * (M - k) > 0.2)) {
      var Y = t[0] - this.center.x, K = t[1] - this.center.y;
      this.radius = Math.sqrt(Y * Y + K * K);
      var l = this.color;
      !l && this.overlay && (l = this.overlay.color), l || (l = "#ff0000"), this.isSelected ? this.selectionColor ? A.strokeStyle = this.selectionColor : A.strokeStyle = gI.increaseBrightness(l, 50) : A.strokeStyle = l, A.lineWidth = this.lineWidth, A.beginPath(), A.arc(this.center.x, this.center.y, this.radius, 0, 2 * Math.PI, !1), C || (this.fillColor && (A.fillStyle = this.fillColor, A.fill()), A.stroke());
    }
  }, B.prototype.isInStroke = function(A, g, C, I) {
    return this.draw(A, g, !0), A.isPointInStroke(C, I);
  }, B.prototype.intersectsBBox = function(A, g, C, I) {
    const E = {
      x: abs(this.center.x - A),
      y: abs(this.center.y - g)
    };
    if (E.x > C / 2 + this.radius || E.y > I / 2 + this.radius)
      return !1;
    if (E.x <= C / 2 || E.y <= I / 2)
      return !0;
    const o = E.x - C / 2, t = E.y - I / 2;
    return o * o + t * t <= this.radius * this.radius;
  }, B;
}(), MQ = function() {
  let B = function(A, g, C, I, E) {
    E = E || {}, this.color = E.color || void 0, this.fillColor = E.fillColor || void 0, this.lineWidth = E.lineWidth || 2, this.id = "ellipse-" + Z.uuidv4(), this.setCenter(A), this.setRadiuses(g, C), this.setRotation(I), this.overlay = null, this.isShowing = !0, this.isSelected = !1, this.selectionColor = "#00ff00";
  };
  return B.prototype.setColor = function(A) {
    this.color != A && (this.color = A, this.overlay && this.overlay.reportChange());
  }, B.prototype.setSelectionColor = function(A) {
    this.selectionColor != A && (this.selectionColor = A, this.overlay && this.overlay.reportChange());
  }, B.prototype.setLineWidth = function(A) {
    this.lineWidth != A && (this.lineWidth = A, this.overlay && this.overlay.reportChange());
  }, B.prototype.setOverlay = function(A) {
    this.overlay = A;
  }, B.prototype.show = function() {
    this.isShowing || (this.isShowing = !0, this.overlay && this.overlay.reportChange());
  }, B.prototype.hide = function() {
    this.isShowing && (this.isShowing = !1, this.overlay && this.overlay.reportChange());
  }, B.prototype.select = function() {
    this.isSelected || (this.isSelected = !0, this.overlay && this.overlay.reportChange());
  }, B.prototype.deselect = function() {
    this.isSelected && (this.isSelected = !1, this.overlay && this.overlay.reportChange());
  }, B.prototype.setCenter = function(A) {
    this.centerRaDec = A, this.overlay && this.overlay.reportChange();
  }, B.prototype.setRotation = function(A) {
    let g = A * Math.PI / 180;
    this.rotation = g, this.overlay && this.overlay.reportChange();
  }, B.prototype.setRadiuses = function(A, g) {
    this.radiusXDegrees = A, this.radiusYDegrees = g, this.overlay && this.overlay.reportChange();
  }, B.prototype.isFootprint = function() {
    return !0;
  }, B.prototype.draw = function(A, g, C) {
    if (!this.isShowing)
      return;
    var I = jA.radecToViewXy(this.centerRaDec[0], this.centerRaDec[1], g);
    if (!I)
      return;
    let E = jA.radecToViewXy(this.centerRaDec[0] + this.radiusXDegrees, this.centerRaDec[1], g), o = jA.radecToViewXy(this.centerRaDec[0], this.centerRaDec[1] + this.radiusYDegrees, g);
    if (!E || !o)
      return;
    var t = E[0] - I[0], s = E[1] - I[1], M = Math.sqrt(t * t + s * s), N = o[0] - I[0], k = o[1] - I[1], U = Math.sqrt(N * N + k * k);
    if (t * k - N * s <= 0)
      return;
    C = C === !0 || !1;
    var Y = this.color;
    !Y && this.overlay && (Y = this.overlay.color), Y || (Y = "#ff0000"), this.isSelected ? this.selectionColor ? A.strokeStyle = this.selectionColor : A.strokeStyle = gI.increaseBrightness(Y, 50) : A.strokeStyle = Y;
    let K = this.centerRaDec, l = [this.centerRaDec[0], this.centerRaDec[1] + 1e-3], d = this.overlay.view.wasm.worldToScreen(K[0], K[1]), q = this.overlay.view.wasm.worldToScreen(l[0], l[1]), W = [q[0] - d[0], q[1] - d[1]], T = Math.sqrt(W[0] * W[0] + W[1] * W[1]);
    W = [W[0] / T, W[1] / T];
    let v = [1, 0], _ = v[0], a = v[1], iA = W[0], IA = W[1], z = Math.atan2(_ * IA - a * iA, _ * iA + a * IA), oA = -this.rotation + z;
    A.beginPath(), A.ellipse(I[0], I[1], M, U, oA, 0, 2 * Math.PI, !1), C || (this.fillColor && (A.fillStyle = this.fillColor, A.fill()), A.stroke());
  }, B.prototype.isInStroke = function(A, g, C, I) {
    return this.draw(A, g, !0), A.isPointInStroke(C, I);
  }, B.prototype.intersectsBBox = function(A, g, C, I) {
  }, B;
}(), uI = function() {
  let B = function(A, g, C, I) {
    this.x1 = A, this.y1 = g, this.x2 = C, this.y2 = I;
  };
  return B.prototype.isInsideView = function(A, g) {
    if (this.x1 >= 0 && this.x1 <= A && this.y1 >= 0 && this.y1 <= g || this.x2 >= 0 && this.x2 <= A && this.y2 >= 0 && this.y2 <= g)
      return !0;
    let C = B.intersectLine(this.x1, this.y1, this.x2, this.y2, 0, 0, 0, g), I = B.intersectLine(this.x1, this.y1, this.x2, this.y2, A, 0, A, g), E = B.intersectLine(this.x1, this.y1, this.x2, this.y2, 0, 0, A, 0), o = B.intersectLine(this.x1, this.y1, this.x2, this.y2, 0, g, A, g);
    return !!(C || I || E || o);
  }, B.prototype.isFootprint = function() {
    return !1;
  }, B.prototype.draw = function(A, g) {
    g = g === !0 || !1, A.beginPath(), A.moveTo(this.x1, this.y1), A.lineTo(this.x2, this.y2), g || A.stroke();
  }, B.intersectLine = function(A, g, C, I, E, o, t, s) {
    let M = ((t - E) * (g - o) - (s - o) * (A - E)) / ((s - o) * (C - A) - (t - E) * (I - g)), N = ((C - A) * (g - o) - (I - g) * (A - E)) / ((s - o) * (C - A) - (t - E) * (I - g));
    return M >= 0 && M <= 1 && N >= 0 && N <= 1;
  }, B.prototype.isInStroke = function(A, g, C, I) {
    return this.draw(A, g, !0), A.isPointInStroke(C, I);
  }, B.prototype.intersectsBBox = function(A, g, C, I) {
  }, B;
}(), eB = function() {
  let B = function(A, g) {
    g = g || {}, this.color = g.color || void 0, this.lineWidth = g.lineWidth || 2, g.closed ? this.closed = g.closed : this.closed = !1, this.id = "polyline-" + Z.uuidv4(), this.radecArray = A, this.overlay = null, this.isShowing = !0, this.isSelected = !1, this.selectionColor = "#00ff00";
  };
  return B.prototype.setOverlay = function(A) {
    this.overlay = A;
  }, B.prototype.show = function() {
    this.isShowing || (this.isShowing = !0, this.overlay && this.overlay.reportChange());
  }, B.prototype.hide = function() {
    this.isShowing && (this.isShowing = !1, this.overlay && this.overlay.reportChange());
  }, B.prototype.select = function() {
    this.isSelected || (this.isSelected = !0, this.overlay && this.overlay.reportChange());
  }, B.prototype.deselect = function() {
    this.isSelected && (this.isSelected = !1, this.overlay && this.overlay.reportChange());
  }, B.prototype.setLineWidth = function(A) {
    this.lineWidth != A && (this.lineWidth = A, this.overlay && this.overlay.reportChange());
  }, B.prototype.setColor = function(A) {
    this.color != A && (this.color = A, this.overlay && this.overlay.reportChange());
  }, B.prototype.setSelectionColor = function(A) {
    this.selectionColor != A && (this.selectionColor = A, this.overlay && this.overlay.reportChange());
  }, B.prototype.isFootprint = function() {
    return this.closed;
  }, B.prototype.draw = function(A, g, C) {
    if (!this.isShowing || !this.radecArray || this.radecArray.length < 2)
      return;
    C = C === !0 || !1;
    var I = this.color;
    !I && this.overlay && (I = this.overlay.color), I || (I = "#ff0000"), this.isSelected ? this.selectionColor ? A.strokeStyle = this.selectionColor : A.strokeStyle = gI.increaseBrightness(I, 50) : A.strokeStyle = I;
    let E = [], o = this.radecArray.length, t = Number.POSITIVE_INFINITY, s = Number.NEGATIVE_INFINITY, M = Number.POSITIVE_INFINITY, N = Number.NEGATIVE_INFINITY;
    for (var k = 0; k < o; k++) {
      var U = jA.radecToViewXy(this.radecArray[k][0], this.radecArray[k][1], g);
      if (!U)
        return;
      E.push({ x: U[0], y: U[1] }), t = Math.min(t, U[0]), M = Math.min(M, U[1]), s = Math.max(s, U[0]), N = Math.max(N, U[1]);
    }
    if (s - t < 1 || N - M < 1)
      return;
    let Y;
    g.projection === uA.SIN ? Y = (q, W) => {
      const T = new uI(q.x, q.y, W.x, W.y);
      T.isInsideView(g.width, g.height) && T.draw(A);
    } : Y = (q, W) => {
      const T = new uI(q.x, q.y, W.x, W.y);
      if (T.isInsideView(g.width, g.height)) {
        const [v, _] = jA.viewXyToClipXy(T.x1, T.y1, g), [a, iA] = jA.viewXyToClipXy(T.x2, T.y2, g);
        (v - a) * (v - a) + (_ - iA) * (_ - iA) < 0.1 && T.draw(A);
      }
    };
    let K = this.closed ? o : o - 1, l = this.closed ? o - 1 : 0, d = this.closed ? 0 : 1;
    A.lineWidth = this.lineWidth, A.beginPath();
    for (var k = 0; k < K; k++)
      Y(E[l], E[d]), l = d, d = d + 1;
    C || A.stroke();
  }, B.prototype.isInStroke = function(A, g, C, I) {
    let E = [];
    for (var o = 0; o < this.radecArray.length; o++) {
      var t = jA.radecToViewXy(this.radecArray[o][0], this.radecArray[o][1], g);
      if (!t)
        return !1;
      E.push({
        x: t[0],
        y: t[1]
      });
    }
    const s = E.length - 1;
    for (var M = 0; M < s; M++)
      if (new uI(E[M].x, E[M].y, E[M + 1].x, E[M + 1].y).draw(A, !0), A.isPointInStroke(C, I))
        return !0;
    return !!(this.closed && (new uI(E[s].x, E[s].y, E[0].x, E[0].y).draw(A, !0), A.isPointInStroke(C, I)));
  }, B.prototype.intersectsBBox = function(A, g, C, I) {
  }, B;
}(), mC = function() {
  return window.requestAnimationFrame || window.webkitRequestAnimationFrame || window.mozRequestAnimationFrame || window.oRequestAnimationFrame || window.msRequestAnimationFrame;
}();
var OC = `#version 300 es
precision lowp float;layout(location=0)in vec2 offset;layout(location=1)in vec2 uv;layout(location=2)in vec3 center;uniform float current_time;uniform mat4 model;uniform mat4 inv_model;uniform vec2 ndc_to_clip;uniform float czf;uniform vec2 kernel_size;out vec2 out_uv;out vec3 out_p;const float PI=3.1415926535897932384626433832795;const mat4 GAL2J2000=mat4(-0.4448296299195045,0.7469822444763707,0.4941094279435681,0.0,-0.1980763734646737,0.4559837762325372,-0.8676661489811610,0.0,-0.873437090247923,-0.4838350155267381,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);const mat4 J20002GAL=mat4(-0.4448296299195045,-0.1980763734646737,-0.873437090247923,0.0,0.7469822444763707,0.4559837762325372,-0.4838350155267381,0.0,0.4941094279435681,-0.8676661489811610,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);vec2 world2clip_orthographic(vec3 p){return vec2(-p.x,p.y);}vec2 world2clip_aitoff(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float theta_by_two=theta*0.5;float alpha=acos(cos(delta)*cos(theta_by_two));float inv_sinc_alpha=1.0;if(alpha>1e-4){inv_sinc_alpha=alpha/sin(alpha);}float x=2.0*inv_sinc_alpha*cos(delta)*sin(theta_by_two);float y=inv_sinc_alpha*sin(delta);return vec2(x/PI,y/PI);}vec2 world2clip_mollweide(vec3 p){int max_iter=10;float delta=asin(p.y);float theta=atan(p.x,p.z);float cst=PI*sin(delta);float phi=delta;float dx=phi+sin(phi)-cst;int k=0;while(abs(dx)>1e-6&&k<max_iter){phi=phi-dx/(1.0+cos(phi));dx=phi+sin(phi)-cst;k=k+1;}phi=phi*0.5;float x=(-theta/PI)*cos(phi);float y=0.5*sin(phi);return vec2(x,y);}vec2 world2clip_mercator(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float x=theta/PI;float y=asinh(tan(delta/PI));return vec2(x,y);}float arc_sinc(float x){if(x>1e-4){return asin(x)/x;}else{float x2=x*x;return 1.0+x2*(1.0+x2*9.0/20.0)/6.0;}}vec2 world2clip_arc(vec3 p){if(p.z>-1.0){float r=length(p.xy);if(p.z>0.0){r=arc_sinc(r);}else{r=acos(p.z)/r;}float x=p.x*r;float y=p.y*r;return vec2(-x/PI,y/PI);}else{return vec2(1.0,0.0);}}vec2 world2clip_gnomonic(vec3 p){if(p.z<=1e-2){return vec2(1.0,0.0);}else{return vec2((-p.x/p.z)/PI,(p.y/p.z)/PI);}}const float TWICE_PI=6.28318530718;const float FOUR_OVER_PI=1.27323954474;const float TRANSITION_Z=0.66666666666;const float TRANSITION_Z_INV=1.5;float one_minus_z_pos(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return 1.0-p.z;}float one_minus_z_neg(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return p.z+1.0;}vec2 xpm1_and_offset(vec2 p){int x_neg=int(p.x<0.0);int y_neg=int(p.y<0.0);int offset=-(y_neg<<2)+1+((x_neg ^ y_neg)<<1);float lon=atan(abs(p.y),abs(p.x));float x02=lon*FOUR_OVER_PI;if(x_neg!=y_neg){return vec2(1.0-x02,float(offset));}else{return vec2(x02-1.0,float(offset));}}vec2 world2clip_healpix(vec3 p){vec2 x_pm1_and_offset=xpm1_and_offset(p.xy);vec2 p_proj=vec2(0.0);if(p.z>TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_pos(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,2.0-sqrt_3_one_min_z);}else if(p.z<-TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_neg(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,-2.0+sqrt_3_one_min_z);}else{p_proj=vec2(atan(p.y,p.x)*FOUR_OVER_PI,p.z*TRANSITION_Z_INV);}return p_proj*vec2(-0.25,0.5);}void main(){vec3 p=vec3(inv_model*vec4(center,1.0f));vec2 center_pos_clip_space=world2clip_aitoff(p);vec2 pos_clip_space=center_pos_clip_space;gl_Position=vec4((pos_clip_space/(ndc_to_clip*czf))+offset*kernel_size,0.f,1.f);out_uv=uv;out_p=p;}`, vC = `#version 300 es
precision lowp float;layout(location=0)in vec2 offset;layout(location=1)in vec2 uv;layout(location=2)in vec3 center;uniform float current_time;uniform mat4 inv_model;uniform vec2 ndc_to_clip;uniform float czf;uniform vec2 kernel_size;out vec2 out_uv;out vec3 out_p;const float PI=3.1415926535897932384626433832795;const mat4 GAL2J2000=mat4(-0.4448296299195045,0.7469822444763707,0.4941094279435681,0.0,-0.1980763734646737,0.4559837762325372,-0.8676661489811610,0.0,-0.873437090247923,-0.4838350155267381,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);const mat4 J20002GAL=mat4(-0.4448296299195045,-0.1980763734646737,-0.873437090247923,0.0,0.7469822444763707,0.4559837762325372,-0.4838350155267381,0.0,0.4941094279435681,-0.8676661489811610,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);vec2 world2clip_orthographic(vec3 p){return vec2(-p.x,p.y);}vec2 world2clip_aitoff(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float theta_by_two=theta*0.5;float alpha=acos(cos(delta)*cos(theta_by_two));float inv_sinc_alpha=1.0;if(alpha>1e-4){inv_sinc_alpha=alpha/sin(alpha);}float x=2.0*inv_sinc_alpha*cos(delta)*sin(theta_by_two);float y=inv_sinc_alpha*sin(delta);return vec2(x/PI,y/PI);}vec2 world2clip_mollweide(vec3 p){int max_iter=10;float delta=asin(p.y);float theta=atan(p.x,p.z);float cst=PI*sin(delta);float phi=delta;float dx=phi+sin(phi)-cst;int k=0;while(abs(dx)>1e-6&&k<max_iter){phi=phi-dx/(1.0+cos(phi));dx=phi+sin(phi)-cst;k=k+1;}phi=phi*0.5;float x=(-theta/PI)*cos(phi);float y=0.5*sin(phi);return vec2(x,y);}vec2 world2clip_mercator(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float x=theta/PI;float y=asinh(tan(delta/PI));return vec2(x,y);}float arc_sinc(float x){if(x>1e-4){return asin(x)/x;}else{float x2=x*x;return 1.0+x2*(1.0+x2*9.0/20.0)/6.0;}}vec2 world2clip_arc(vec3 p){if(p.z>-1.0){float r=length(p.xy);if(p.z>0.0){r=arc_sinc(r);}else{r=acos(p.z)/r;}float x=p.x*r;float y=p.y*r;return vec2(-x/PI,y/PI);}else{return vec2(1.0,0.0);}}vec2 world2clip_gnomonic(vec3 p){if(p.z<=1e-2){return vec2(1.0,0.0);}else{return vec2((-p.x/p.z)/PI,(p.y/p.z)/PI);}}const float TWICE_PI=6.28318530718;const float FOUR_OVER_PI=1.27323954474;const float TRANSITION_Z=0.66666666666;const float TRANSITION_Z_INV=1.5;float one_minus_z_pos(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return 1.0-p.z;}float one_minus_z_neg(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return p.z+1.0;}vec2 xpm1_and_offset(vec2 p){int x_neg=int(p.x<0.0);int y_neg=int(p.y<0.0);int offset=-(y_neg<<2)+1+((x_neg ^ y_neg)<<1);float lon=atan(abs(p.y),abs(p.x));float x02=lon*FOUR_OVER_PI;if(x_neg!=y_neg){return vec2(1.0-x02,float(offset));}else{return vec2(x02-1.0,float(offset));}}vec2 world2clip_healpix(vec3 p){vec2 x_pm1_and_offset=xpm1_and_offset(p.xy);vec2 p_proj=vec2(0.0);if(p.z>TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_pos(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,2.0-sqrt_3_one_min_z);}else if(p.z<-TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_neg(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,-2.0+sqrt_3_one_min_z);}else{p_proj=vec2(atan(p.y,p.x)*FOUR_OVER_PI,p.z*TRANSITION_Z_INV);}return p_proj*vec2(-0.25,0.5);}void main(){vec3 p=vec3(inv_model*vec4(center,1.0f));vec2 center_pos_clip_space=world2clip_mercator(p);vec2 pos_clip_space=center_pos_clip_space;gl_Position=vec4((pos_clip_space/(ndc_to_clip*czf))+offset*kernel_size,0.f,1.f);out_uv=uv;out_p=p;}`, bC = `#version 300 es
precision lowp float;layout(location=0)in vec2 offset;layout(location=1)in vec2 uv;layout(location=2)in vec3 center;uniform float current_time;uniform mat4 inv_model;uniform vec2 ndc_to_clip;uniform float czf;uniform vec2 kernel_size;out vec2 out_uv;out vec3 out_p;const float PI=3.1415926535897932384626433832795;const mat4 GAL2J2000=mat4(-0.4448296299195045,0.7469822444763707,0.4941094279435681,0.0,-0.1980763734646737,0.4559837762325372,-0.8676661489811610,0.0,-0.873437090247923,-0.4838350155267381,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);const mat4 J20002GAL=mat4(-0.4448296299195045,-0.1980763734646737,-0.873437090247923,0.0,0.7469822444763707,0.4559837762325372,-0.4838350155267381,0.0,0.4941094279435681,-0.8676661489811610,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);vec2 world2clip_orthographic(vec3 p){return vec2(-p.x,p.y);}vec2 world2clip_aitoff(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float theta_by_two=theta*0.5;float alpha=acos(cos(delta)*cos(theta_by_two));float inv_sinc_alpha=1.0;if(alpha>1e-4){inv_sinc_alpha=alpha/sin(alpha);}float x=2.0*inv_sinc_alpha*cos(delta)*sin(theta_by_two);float y=inv_sinc_alpha*sin(delta);return vec2(x/PI,y/PI);}vec2 world2clip_mollweide(vec3 p){int max_iter=10;float delta=asin(p.y);float theta=atan(p.x,p.z);float cst=PI*sin(delta);float phi=delta;float dx=phi+sin(phi)-cst;int k=0;while(abs(dx)>1e-6&&k<max_iter){phi=phi-dx/(1.0+cos(phi));dx=phi+sin(phi)-cst;k=k+1;}phi=phi*0.5;float x=(-theta/PI)*cos(phi);float y=0.5*sin(phi);return vec2(x,y);}vec2 world2clip_mercator(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float x=theta/PI;float y=asinh(tan(delta/PI));return vec2(x,y);}float arc_sinc(float x){if(x>1e-4){return asin(x)/x;}else{float x2=x*x;return 1.0+x2*(1.0+x2*9.0/20.0)/6.0;}}vec2 world2clip_arc(vec3 p){if(p.z>-1.0){float r=length(p.xy);if(p.z>0.0){r=arc_sinc(r);}else{r=acos(p.z)/r;}float x=p.x*r;float y=p.y*r;return vec2(-x/PI,y/PI);}else{return vec2(1.0,0.0);}}vec2 world2clip_gnomonic(vec3 p){if(p.z<=1e-2){return vec2(1.0,0.0);}else{return vec2((-p.x/p.z)/PI,(p.y/p.z)/PI);}}const float TWICE_PI=6.28318530718;const float FOUR_OVER_PI=1.27323954474;const float TRANSITION_Z=0.66666666666;const float TRANSITION_Z_INV=1.5;float one_minus_z_pos(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return 1.0-p.z;}float one_minus_z_neg(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return p.z+1.0;}vec2 xpm1_and_offset(vec2 p){int x_neg=int(p.x<0.0);int y_neg=int(p.y<0.0);int offset=-(y_neg<<2)+1+((x_neg ^ y_neg)<<1);float lon=atan(abs(p.y),abs(p.x));float x02=lon*FOUR_OVER_PI;if(x_neg!=y_neg){return vec2(1.0-x02,float(offset));}else{return vec2(x02-1.0,float(offset));}}vec2 world2clip_healpix(vec3 p){vec2 x_pm1_and_offset=xpm1_and_offset(p.xy);vec2 p_proj=vec2(0.0);if(p.z>TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_pos(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,2.0-sqrt_3_one_min_z);}else if(p.z<-TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_neg(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,-2.0+sqrt_3_one_min_z);}else{p_proj=vec2(atan(p.y,p.x)*FOUR_OVER_PI,p.z*TRANSITION_Z_INV);}return p_proj*vec2(-0.25,0.5);}void main(){vec3 p=vec3(inv_model*vec4(center,1.0f));vec2 center_pos_clip_space=world2clip_arc(p);vec2 pos_clip_space=center_pos_clip_space;gl_Position=vec4((pos_clip_space/(ndc_to_clip*czf))+offset*kernel_size,0.f,1.f);out_uv=uv;out_p=p;}`, ZC = `#version 300 es
precision lowp float;layout(location=0)in vec2 offset;layout(location=1)in vec2 uv;layout(location=2)in vec3 center;uniform float current_time;uniform mat4 inv_model;uniform vec2 ndc_to_clip;uniform float czf;uniform vec2 kernel_size;out vec2 out_uv;out vec3 out_p;const float PI=3.1415926535897932384626433832795;const mat4 GAL2J2000=mat4(-0.4448296299195045,0.7469822444763707,0.4941094279435681,0.0,-0.1980763734646737,0.4559837762325372,-0.8676661489811610,0.0,-0.873437090247923,-0.4838350155267381,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);const mat4 J20002GAL=mat4(-0.4448296299195045,-0.1980763734646737,-0.873437090247923,0.0,0.7469822444763707,0.4559837762325372,-0.4838350155267381,0.0,0.4941094279435681,-0.8676661489811610,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);vec2 world2clip_orthographic(vec3 p){return vec2(-p.x,p.y);}vec2 world2clip_aitoff(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float theta_by_two=theta*0.5;float alpha=acos(cos(delta)*cos(theta_by_two));float inv_sinc_alpha=1.0;if(alpha>1e-4){inv_sinc_alpha=alpha/sin(alpha);}float x=2.0*inv_sinc_alpha*cos(delta)*sin(theta_by_two);float y=inv_sinc_alpha*sin(delta);return vec2(x/PI,y/PI);}vec2 world2clip_mollweide(vec3 p){int max_iter=10;float delta=asin(p.y);float theta=atan(p.x,p.z);float cst=PI*sin(delta);float phi=delta;float dx=phi+sin(phi)-cst;int k=0;while(abs(dx)>1e-6&&k<max_iter){phi=phi-dx/(1.0+cos(phi));dx=phi+sin(phi)-cst;k=k+1;}phi=phi*0.5;float x=(-theta/PI)*cos(phi);float y=0.5*sin(phi);return vec2(x,y);}vec2 world2clip_mercator(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float x=theta/PI;float y=asinh(tan(delta/PI));return vec2(x,y);}float arc_sinc(float x){if(x>1e-4){return asin(x)/x;}else{float x2=x*x;return 1.0+x2*(1.0+x2*9.0/20.0)/6.0;}}vec2 world2clip_arc(vec3 p){if(p.z>-1.0){float r=length(p.xy);if(p.z>0.0){r=arc_sinc(r);}else{r=acos(p.z)/r;}float x=p.x*r;float y=p.y*r;return vec2(-x/PI,y/PI);}else{return vec2(1.0,0.0);}}vec2 world2clip_gnomonic(vec3 p){if(p.z<=1e-2){return vec2(1.0,0.0);}else{return vec2((-p.x/p.z)/PI,(p.y/p.z)/PI);}}const float TWICE_PI=6.28318530718;const float FOUR_OVER_PI=1.27323954474;const float TRANSITION_Z=0.66666666666;const float TRANSITION_Z_INV=1.5;float one_minus_z_pos(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return 1.0-p.z;}float one_minus_z_neg(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return p.z+1.0;}vec2 xpm1_and_offset(vec2 p){int x_neg=int(p.x<0.0);int y_neg=int(p.y<0.0);int offset=-(y_neg<<2)+1+((x_neg ^ y_neg)<<1);float lon=atan(abs(p.y),abs(p.x));float x02=lon*FOUR_OVER_PI;if(x_neg!=y_neg){return vec2(1.0-x02,float(offset));}else{return vec2(x02-1.0,float(offset));}}vec2 world2clip_healpix(vec3 p){vec2 x_pm1_and_offset=xpm1_and_offset(p.xy);vec2 p_proj=vec2(0.0);if(p.z>TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_pos(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,2.0-sqrt_3_one_min_z);}else if(p.z<-TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_neg(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,-2.0+sqrt_3_one_min_z);}else{p_proj=vec2(atan(p.y,p.x)*FOUR_OVER_PI,p.z*TRANSITION_Z_INV);}return p_proj*vec2(-0.25,0.5);}void main(){vec3 p=vec3(inv_model*vec4(center,1.0f));vec2 center_pos_clip_space=world2clip_gnomonic(p);vec2 pos_clip_space=center_pos_clip_space;gl_Position=vec4((pos_clip_space/(ndc_to_clip*czf))+offset*kernel_size,0.f,1.f);out_uv=uv;out_p=p;}`, TC = `#version 300 es
precision lowp float;layout(location=0)in vec2 offset;layout(location=1)in vec2 uv;layout(location=2)in vec3 center;uniform float current_time;uniform mat4 inv_model;uniform vec2 ndc_to_clip;uniform float czf;uniform vec2 kernel_size;out vec2 out_uv;out vec3 out_p;const float PI=3.1415926535897932384626433832795;const mat4 GAL2J2000=mat4(-0.4448296299195045,0.7469822444763707,0.4941094279435681,0.0,-0.1980763734646737,0.4559837762325372,-0.8676661489811610,0.0,-0.873437090247923,-0.4838350155267381,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);const mat4 J20002GAL=mat4(-0.4448296299195045,-0.1980763734646737,-0.873437090247923,0.0,0.7469822444763707,0.4559837762325372,-0.4838350155267381,0.0,0.4941094279435681,-0.8676661489811610,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);vec2 world2clip_orthographic(vec3 p){return vec2(-p.x,p.y);}vec2 world2clip_aitoff(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float theta_by_two=theta*0.5;float alpha=acos(cos(delta)*cos(theta_by_two));float inv_sinc_alpha=1.0;if(alpha>1e-4){inv_sinc_alpha=alpha/sin(alpha);}float x=2.0*inv_sinc_alpha*cos(delta)*sin(theta_by_two);float y=inv_sinc_alpha*sin(delta);return vec2(x/PI,y/PI);}vec2 world2clip_mollweide(vec3 p){int max_iter=10;float delta=asin(p.y);float theta=atan(p.x,p.z);float cst=PI*sin(delta);float phi=delta;float dx=phi+sin(phi)-cst;int k=0;while(abs(dx)>1e-6&&k<max_iter){phi=phi-dx/(1.0+cos(phi));dx=phi+sin(phi)-cst;k=k+1;}phi=phi*0.5;float x=(-theta/PI)*cos(phi);float y=0.5*sin(phi);return vec2(x,y);}vec2 world2clip_mercator(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float x=theta/PI;float y=asinh(tan(delta/PI));return vec2(x,y);}float arc_sinc(float x){if(x>1e-4){return asin(x)/x;}else{float x2=x*x;return 1.0+x2*(1.0+x2*9.0/20.0)/6.0;}}vec2 world2clip_arc(vec3 p){if(p.z>-1.0){float r=length(p.xy);if(p.z>0.0){r=arc_sinc(r);}else{r=acos(p.z)/r;}float x=p.x*r;float y=p.y*r;return vec2(-x/PI,y/PI);}else{return vec2(1.0,0.0);}}vec2 world2clip_gnomonic(vec3 p){if(p.z<=1e-2){return vec2(1.0,0.0);}else{return vec2((-p.x/p.z)/PI,(p.y/p.z)/PI);}}const float TWICE_PI=6.28318530718;const float FOUR_OVER_PI=1.27323954474;const float TRANSITION_Z=0.66666666666;const float TRANSITION_Z_INV=1.5;float one_minus_z_pos(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return 1.0-p.z;}float one_minus_z_neg(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return p.z+1.0;}vec2 xpm1_and_offset(vec2 p){int x_neg=int(p.x<0.0);int y_neg=int(p.y<0.0);int offset=-(y_neg<<2)+1+((x_neg ^ y_neg)<<1);float lon=atan(abs(p.y),abs(p.x));float x02=lon*FOUR_OVER_PI;if(x_neg!=y_neg){return vec2(1.0-x02,float(offset));}else{return vec2(x02-1.0,float(offset));}}vec2 world2clip_healpix(vec3 p){vec2 x_pm1_and_offset=xpm1_and_offset(p.xy);vec2 p_proj=vec2(0.0);if(p.z>TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_pos(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,2.0-sqrt_3_one_min_z);}else if(p.z<-TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_neg(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,-2.0+sqrt_3_one_min_z);}else{p_proj=vec2(atan(p.y,p.x)*FOUR_OVER_PI,p.z*TRANSITION_Z_INV);}return p_proj*vec2(-0.25,0.5);}void main(){vec3 p=vec3(inv_model*vec4(center,1.0f));vec2 center_pos_clip_space=world2clip_mollweide(p);vec2 pos_clip_space=center_pos_clip_space;gl_Position=vec4((pos_clip_space/(ndc_to_clip*czf))+offset*kernel_size,0.f,1.f);out_uv=uv;out_p=p;}`, WC = `#version 300 es
precision lowp float;layout(location=0)in vec2 offset;layout(location=1)in vec2 uv;layout(location=2)in vec3 center;uniform float current_time;uniform mat4 model;uniform mat4 inv_model;uniform vec2 ndc_to_clip;uniform float czf;uniform vec2 kernel_size;out vec2 out_uv;out vec3 out_p;const float PI=3.1415926535897932384626433832795;const mat4 GAL2J2000=mat4(-0.4448296299195045,0.7469822444763707,0.4941094279435681,0.0,-0.1980763734646737,0.4559837762325372,-0.8676661489811610,0.0,-0.873437090247923,-0.4838350155267381,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);const mat4 J20002GAL=mat4(-0.4448296299195045,-0.1980763734646737,-0.873437090247923,0.0,0.7469822444763707,0.4559837762325372,-0.4838350155267381,0.0,0.4941094279435681,-0.8676661489811610,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);vec2 world2clip_orthographic(vec3 p){return vec2(-p.x,p.y);}vec2 world2clip_aitoff(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float theta_by_two=theta*0.5;float alpha=acos(cos(delta)*cos(theta_by_two));float inv_sinc_alpha=1.0;if(alpha>1e-4){inv_sinc_alpha=alpha/sin(alpha);}float x=2.0*inv_sinc_alpha*cos(delta)*sin(theta_by_two);float y=inv_sinc_alpha*sin(delta);return vec2(x/PI,y/PI);}vec2 world2clip_mollweide(vec3 p){int max_iter=10;float delta=asin(p.y);float theta=atan(p.x,p.z);float cst=PI*sin(delta);float phi=delta;float dx=phi+sin(phi)-cst;int k=0;while(abs(dx)>1e-6&&k<max_iter){phi=phi-dx/(1.0+cos(phi));dx=phi+sin(phi)-cst;k=k+1;}phi=phi*0.5;float x=(-theta/PI)*cos(phi);float y=0.5*sin(phi);return vec2(x,y);}vec2 world2clip_mercator(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float x=theta/PI;float y=asinh(tan(delta/PI));return vec2(x,y);}float arc_sinc(float x){if(x>1e-4){return asin(x)/x;}else{float x2=x*x;return 1.0+x2*(1.0+x2*9.0/20.0)/6.0;}}vec2 world2clip_arc(vec3 p){if(p.z>-1.0){float r=length(p.xy);if(p.z>0.0){r=arc_sinc(r);}else{r=acos(p.z)/r;}float x=p.x*r;float y=p.y*r;return vec2(-x/PI,y/PI);}else{return vec2(1.0,0.0);}}vec2 world2clip_gnomonic(vec3 p){if(p.z<=1e-2){return vec2(1.0,0.0);}else{return vec2((-p.x/p.z)/PI,(p.y/p.z)/PI);}}const float TWICE_PI=6.28318530718;const float FOUR_OVER_PI=1.27323954474;const float TRANSITION_Z=0.66666666666;const float TRANSITION_Z_INV=1.5;float one_minus_z_pos(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return 1.0-p.z;}float one_minus_z_neg(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return p.z+1.0;}vec2 xpm1_and_offset(vec2 p){int x_neg=int(p.x<0.0);int y_neg=int(p.y<0.0);int offset=-(y_neg<<2)+1+((x_neg ^ y_neg)<<1);float lon=atan(abs(p.y),abs(p.x));float x02=lon*FOUR_OVER_PI;if(x_neg!=y_neg){return vec2(1.0-x02,float(offset));}else{return vec2(x02-1.0,float(offset));}}vec2 world2clip_healpix(vec3 p){vec2 x_pm1_and_offset=xpm1_and_offset(p.xy);vec2 p_proj=vec2(0.0);if(p.z>TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_pos(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,2.0-sqrt_3_one_min_z);}else if(p.z<-TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_neg(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,-2.0+sqrt_3_one_min_z);}else{p_proj=vec2(atan(p.y,p.x)*FOUR_OVER_PI,p.z*TRANSITION_Z_INV);}return p_proj*vec2(-0.25,0.5);}void main(){vec3 p=vec3(inv_model*vec4(center,1.0f));vec2 center_pos_clip_space=world2clip_healpix(p);vec2 pos_clip_space=center_pos_clip_space;gl_Position=vec4((pos_clip_space/(ndc_to_clip*czf))+offset*kernel_size,0.f,1.f);out_uv=uv;out_p=p;}`, VC = `#version 300 es
precision lowp float;layout(location=0)in vec2 offset;layout(location=1)in vec2 uv;layout(location=2)in vec3 center;uniform float current_time;uniform mat4 inv_model;uniform vec2 ndc_to_clip;uniform float czf;uniform vec2 kernel_size;out vec2 out_uv;out vec3 out_p;const float PI=3.1415926535897932384626433832795;const mat4 GAL2J2000=mat4(-0.4448296299195045,0.7469822444763707,0.4941094279435681,0.0,-0.1980763734646737,0.4559837762325372,-0.8676661489811610,0.0,-0.873437090247923,-0.4838350155267381,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);const mat4 J20002GAL=mat4(-0.4448296299195045,-0.1980763734646737,-0.873437090247923,0.0,0.7469822444763707,0.4559837762325372,-0.4838350155267381,0.0,0.4941094279435681,-0.8676661489811610,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);vec2 world2clip_orthographic(vec3 p){return vec2(-p.x,p.y);}vec2 world2clip_aitoff(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float theta_by_two=theta*0.5;float alpha=acos(cos(delta)*cos(theta_by_two));float inv_sinc_alpha=1.0;if(alpha>1e-4){inv_sinc_alpha=alpha/sin(alpha);}float x=2.0*inv_sinc_alpha*cos(delta)*sin(theta_by_two);float y=inv_sinc_alpha*sin(delta);return vec2(x/PI,y/PI);}vec2 world2clip_mollweide(vec3 p){int max_iter=10;float delta=asin(p.y);float theta=atan(p.x,p.z);float cst=PI*sin(delta);float phi=delta;float dx=phi+sin(phi)-cst;int k=0;while(abs(dx)>1e-6&&k<max_iter){phi=phi-dx/(1.0+cos(phi));dx=phi+sin(phi)-cst;k=k+1;}phi=phi*0.5;float x=(-theta/PI)*cos(phi);float y=0.5*sin(phi);return vec2(x,y);}vec2 world2clip_mercator(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float x=theta/PI;float y=asinh(tan(delta/PI));return vec2(x,y);}float arc_sinc(float x){if(x>1e-4){return asin(x)/x;}else{float x2=x*x;return 1.0+x2*(1.0+x2*9.0/20.0)/6.0;}}vec2 world2clip_arc(vec3 p){if(p.z>-1.0){float r=length(p.xy);if(p.z>0.0){r=arc_sinc(r);}else{r=acos(p.z)/r;}float x=p.x*r;float y=p.y*r;return vec2(-x/PI,y/PI);}else{return vec2(1.0,0.0);}}vec2 world2clip_gnomonic(vec3 p){if(p.z<=1e-2){return vec2(1.0,0.0);}else{return vec2((-p.x/p.z)/PI,(p.y/p.z)/PI);}}const float TWICE_PI=6.28318530718;const float FOUR_OVER_PI=1.27323954474;const float TRANSITION_Z=0.66666666666;const float TRANSITION_Z_INV=1.5;float one_minus_z_pos(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return 1.0-p.z;}float one_minus_z_neg(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return p.z+1.0;}vec2 xpm1_and_offset(vec2 p){int x_neg=int(p.x<0.0);int y_neg=int(p.y<0.0);int offset=-(y_neg<<2)+1+((x_neg ^ y_neg)<<1);float lon=atan(abs(p.y),abs(p.x));float x02=lon*FOUR_OVER_PI;if(x_neg!=y_neg){return vec2(1.0-x02,float(offset));}else{return vec2(x02-1.0,float(offset));}}vec2 world2clip_healpix(vec3 p){vec2 x_pm1_and_offset=xpm1_and_offset(p.xy);vec2 p_proj=vec2(0.0);if(p.z>TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_pos(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,2.0-sqrt_3_one_min_z);}else if(p.z<-TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_neg(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,-2.0+sqrt_3_one_min_z);}else{p_proj=vec2(atan(p.y,p.x)*FOUR_OVER_PI,p.z*TRANSITION_Z_INV);}return p_proj*vec2(-0.25,0.5);}void main(){vec3 p=vec3(inv_model*vec4(center,1.0f));vec2 center_pos_clip_space=world2clip_orthographic(p);vec2 pos_clip_space=center_pos_clip_space;gl_Position=vec4((pos_clip_space/(ndc_to_clip*czf))+offset*kernel_size,0.f,1.f);out_uv=uv;out_p=p;}`, jC = `#version 300 es
precision lowp float;in vec2 out_uv;in vec3 out_p;out vec4 color;uniform sampler2D kernel_texture;uniform float max_density;uniform float fov;uniform float strength;void main(){if(out_p.z<0.f){discard;}color=texture(kernel_texture,out_uv)/max(log2(fov*100.0),1.0);color.r*=strength;}`, _C = `#version 300 es
precision lowp float;in vec2 out_uv;in vec3 out_p;out vec4 color;uniform sampler2D kernel_texture;uniform float max_density;uniform float fov;uniform float strength;void main(){color=texture(kernel_texture,out_uv)/max(log2(fov*100.0),1.0);color.r*=strength;}`, PC = `#version 300 es
precision lowp float;precision lowp sampler2D;layout(location=0)in vec2 position;layout(location=1)in vec2 uv;out vec2 out_uv;void main(){gl_Position=vec4(position,0.f,1.f);out_uv=uv;}`, XC = `#version 300 es
precision lowp float;precision lowp sampler2D;in vec2 out_uv;out vec4 color;uniform sampler2D texture_fbo;uniform float alpha;uniform sampler2D colormaps;uniform float num_colormaps;uniform float colormap_id;vec4 colormap_f(float x){float id=(colormap_id+0.5)/num_colormaps;return texture(colormaps,vec2(x,id));}void main(){float opacity=texture(texture_fbo,out_uv).r;float o=smoothstep(0.0,0.1,opacity);color=colormap_f(opacity);color.a=o*alpha;}`, zC = `#version 300 es
precision lowp float;layout(location=0)in vec2 ndc_pos;void main(){gl_Position=vec4(ndc_pos,0.0,1.0);}`, $C = `#version 300 es
precision lowp float;out vec4 frag_color;uniform vec3 color;uniform float opacity;const float PI=3.141592653589793f;void main(){frag_color=vec4(color,opacity);}`, AE = `#version 300 es
precision highp float;precision mediump int;layout(location=0)in vec2 pos_clip_space;out vec2 out_clip_pos;out vec3 frag_pos;uniform vec2 ndc_to_clip;uniform float czf;uniform mat4 model;uniform sampler2D position_tex;const float PI=3.1415926535897932384626433832795;const mat4 GAL2J2000=mat4(-0.4448296299195045,0.7469822444763707,0.4941094279435681,0.0,-0.1980763734646737,0.4559837762325372,-0.8676661489811610,0.0,-0.873437090247923,-0.4838350155267381,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);const mat4 J20002GAL=mat4(-0.4448296299195045,-0.1980763734646737,-0.873437090247923,0.0,0.7469822444763707,0.4559837762325372,-0.4838350155267381,0.0,0.4941094279435681,-0.8676661489811610,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);vec2 world2clip_orthographic(vec3 p){return vec2(-p.x,p.y);}vec2 world2clip_aitoff(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float theta_by_two=theta*0.5;float alpha=acos(cos(delta)*cos(theta_by_two));float inv_sinc_alpha=1.0;if(alpha>1e-4){inv_sinc_alpha=alpha/sin(alpha);}float x=2.0*inv_sinc_alpha*cos(delta)*sin(theta_by_two);float y=inv_sinc_alpha*sin(delta);return vec2(x/PI,y/PI);}vec2 world2clip_mollweide(vec3 p){int max_iter=10;float delta=asin(p.y);float theta=atan(p.x,p.z);float cst=PI*sin(delta);float phi=delta;float dx=phi+sin(phi)-cst;int k=0;while(abs(dx)>1e-6&&k<max_iter){phi=phi-dx/(1.0+cos(phi));dx=phi+sin(phi)-cst;k=k+1;}phi=phi*0.5;float x=(-theta/PI)*cos(phi);float y=0.5*sin(phi);return vec2(x,y);}vec2 world2clip_mercator(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float x=theta/PI;float y=asinh(tan(delta/PI));return vec2(x,y);}float arc_sinc(float x){if(x>1e-4){return asin(x)/x;}else{float x2=x*x;return 1.0+x2*(1.0+x2*9.0/20.0)/6.0;}}vec2 world2clip_arc(vec3 p){if(p.z>-1.0){float r=length(p.xy);if(p.z>0.0){r=arc_sinc(r);}else{r=acos(p.z)/r;}float x=p.x*r;float y=p.y*r;return vec2(-x/PI,y/PI);}else{return vec2(1.0,0.0);}}vec2 world2clip_gnomonic(vec3 p){if(p.z<=1e-2){return vec2(1.0,0.0);}else{return vec2((-p.x/p.z)/PI,(p.y/p.z)/PI);}}const float TWICE_PI=6.28318530718;const float FOUR_OVER_PI=1.27323954474;const float TRANSITION_Z=0.66666666666;const float TRANSITION_Z_INV=1.5;float one_minus_z_pos(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return 1.0-p.z;}float one_minus_z_neg(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return p.z+1.0;}vec2 xpm1_and_offset(vec2 p){int x_neg=int(p.x<0.0);int y_neg=int(p.y<0.0);int offset=-(y_neg<<2)+1+((x_neg ^ y_neg)<<1);float lon=atan(abs(p.y),abs(p.x));float x02=lon*FOUR_OVER_PI;if(x_neg!=y_neg){return vec2(1.0-x02,float(offset));}else{return vec2(x02-1.0,float(offset));}}vec2 world2clip_healpix(vec3 p){vec2 x_pm1_and_offset=xpm1_and_offset(p.xy);vec2 p_proj=vec2(0.0);if(p.z>TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_pos(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,2.0-sqrt_3_one_min_z);}else if(p.z<-TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_neg(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,-2.0+sqrt_3_one_min_z);}else{p_proj=vec2(atan(p.y,p.x)*FOUR_OVER_PI,p.z*TRANSITION_Z_INV);}return p_proj*vec2(-0.25,0.5);}void main(){vec2 uv=pos_clip_space*0.5+0.5;vec3 world_pos=texture(position_tex,uv).rgb;frag_pos=vec3(model*vec4(world_pos,1.0));gl_Position=vec4(pos_clip_space/(ndc_to_clip*czf),0.0,1.0);out_clip_pos=pos_clip_space;}`, gE = `#version 300 es
precision highp float;precision highp sampler2D;precision highp usampler2D;precision highp isampler2D;precision mediump int;in vec2 out_clip_pos;in vec3 frag_pos;out vec4 out_frag_color;struct Tile{int uniq;int texture_idx;float start_time;float empty;};uniform Tile textures_tiles[12];uniform int num_tiles;uniform sampler2D tex1;uniform sampler2D tex2;uniform sampler2D tex3;uniform int num_tex;uniform float scale;uniform float offset;uniform float blank;uniform float min_value;uniform float max_value;uniform int H;uniform float size_tile_uv;uniform int tex_storing_fits;uniform sampler2D colormaps;uniform float num_colormaps;uniform float colormap_id;vec4 colormap_f(float x){float id=(colormap_id+0.5)/num_colormaps;return texture(colormaps,vec2(x,id));}float linear_f(float x,float min_value,float max_value){return clamp((x-min_value)/(max_value-min_value),0.0,1.0);}float sqrt_f(float x,float min_value,float max_value){float a=linear_f(x,min_value,max_value);return sqrt(a);}float log_f(float x,float min_value,float max_value){float y=linear_f(x,min_value,max_value);float a=1000.0;return log(a*y+1.0)/log(a);}float asinh_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return asinh(10.0*d)/3.0;}float pow2_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return d*d;}float transfer_func(int H,float x,float min_value,float max_value){if(H==0){return linear_f(x,min_value,max_value);}else if(H==1){return sqrt_f(x,min_value,max_value);}else if(H==2){return log_f(x,min_value,max_value);}else if(H==3){return asinh_f(x,min_value,max_value);}else{return pow2_f(x,min_value,max_value);}}uniform float k_gamma;uniform float k_saturation;uniform float k_contrast;uniform float k_brightness;uniform float k_exposure;vec4 apply_gamma(vec4 ic,float g){float new_r=pow(ic.r,g);float new_g=pow(ic.g,g);float new_b=pow(ic.b,g);return vec4(new_r,new_g,new_b,ic.a);}vec4 apply_saturation(vec4 color,float value){const vec3 luminosity_factor=vec3(0.2126,0.7152,0.0722);vec3 grayscale=vec3(dot(color.rgb,luminosity_factor));return vec4(mix(grayscale,color.rgb,1.0+value),color.a);}vec4 apply_contrast(vec4 color,float value){return vec4(0.5+(1.0+value)*(color.rgb-0.5),color.a);}vec4 apply_brightness(vec4 color,float value){return vec4(color.rgb+value,color.a);}vec4 apply_exposure(vec4 color,float value){return vec4((1.0+value)*color.rgb,color.a);}vec4 apply_tonal(vec4 color){return apply_gamma(apply_saturation(apply_contrast(apply_brightness(color,k_brightness),k_contrast),k_saturation),k_gamma);}vec3 rgb2hsv(vec3 c){vec4 K=vec4(0.0,-1.0/3.0,2.0/3.0,-1.0);vec4 p=mix(vec4(c.bg,K.wz),vec4(c.gb,K.xy),step(c.b,c.g));vec4 q=mix(vec4(p.xyw,c.r),vec4(c.r,p.yzx),step(p.x,c.r));float d=q.x-min(q.w,q.y);float e=1.0e-10;return vec3(abs(q.z+(q.w-q.y)/(6.0*d+e)),d/(q.x+e),q.x);}vec3 hsv2rgb(vec3 c){vec4 K=vec4(1.0,2.0/3.0,1.0/3.0,3.0);vec3 p=abs(fract(c.xxx+K.xyz)*6.0-K.www);return c.z*mix(K.xxx,clamp(p-K.xxx,0.0,1.0),c.y);}vec4 get_pixels(vec3 uv){int idx_texture=int(uv.z);if(idx_texture==0){return texture(tex1,uv.xy);}else if(idx_texture==1){return texture(tex2,uv.xy);}else if(idx_texture==2){return texture(tex3,uv.xy);}else{return vec4(0.0,1.0,0.0,1.0);}}vec3 reverse_uv(vec3 uv){uv.y=size_tile_uv+2.0*size_tile_uv*floor(uv.y/size_tile_uv)-uv.y;return uv;}uniform float reversed;vec4 get_color_from_texture(vec3 UV){vec4 color=get_pixels(UV);color.r=transfer_func(H,color.r,min_value,max_value);color.g=transfer_func(H,color.g,min_value,max_value);color.b=transfer_func(H,color.b,min_value,max_value);color.rgb=mix(color.rgb,1.0-color.rgb,reversed);return apply_tonal(color);}vec4 apply_colormap_to_grayscale(float x,float a){float alpha=x*scale+offset;alpha=transfer_func(H,alpha,min_value,max_value);alpha=mix(alpha,1.0-alpha,reversed);vec4 new_color=mix(colormap_f(alpha)*a,vec4(0.0),float(x==blank||isnan(x)));return apply_tonal(new_color);}vec4 get_colormap_from_grayscale_texture(vec3 UV){vec3 uv=mix(UV,reverse_uv(UV),float(tex_storing_fits==1));vec4 color=get_pixels(uv);return apply_colormap_to_grayscale(color.r,color.a);}const float TWICE_PI=6.28318530718;const float PI=3.141592653589793;const float FOUR_OVER_PI=1.27323954474;const float TRANSITION_Z=0.66666666666;const float TRANSITION_Z_INV=1.5;int quarter(vec2 p){int x_neg=int(p.x<0.0);int y_neg=int(p.y<0.0);int q=(x_neg+y_neg)|(y_neg<<1);return q;}float xpm1(vec2 p){bool x_neg=(p.x<0.0);bool y_neg=(p.y<0.0);float lon=atan(abs(p.y),abs(p.x));float x02=lon*FOUR_OVER_PI;if(x_neg!=y_neg){return 1.0-x02;}else{return x02-1.0;}}float one_minus_z_pos(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return 1.0f-p.z;}float one_minus_z_neg(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1f){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return p.z+1.0;}int ij2z(int i,int j){int i4=i|(j<<2);int j4=(i4 ^(i4>>1))&0x22222222;int i5=i4 ^ j4 ^(j4<<1);return i5;}struct HashDxDy{int idx;float dx;float dy;};uniform sampler2D ang2pixd;HashDxDy hash_with_dxdy2(vec2 radec){vec2 aa=vec2(radec.x/TWICE_PI+1.0,(radec.y/PI)+0.5);vec3 v=texture(ang2pixd,aa).rgb;return HashDxDy(int(v.x*255.0),v.y,v.z);}HashDxDy hash_with_dxdy(int depth,vec3 p){int nside=1<<depth;float half_nside=float(nside)*0.5;float x_pm1=xpm1(p.xy);int q=quarter(p.xy);int d0h=0;vec2 p_proj=vec2(0.0);if(p.z>TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_pos(p));p_proj=vec2(x_pm1*sqrt_3_one_min_z,2.0-sqrt_3_one_min_z);d0h=q;}else if(p.z<-TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_neg(p));p_proj=vec2(x_pm1*sqrt_3_one_min_z,sqrt_3_one_min_z);d0h=q+8;}else{float y_pm1=p.z*TRANSITION_Z_INV;int q01=int(x_pm1>y_pm1);int q12=int(x_pm1>=-y_pm1);int q03=1-q12;int q1=q01&q12;p_proj=vec2(x_pm1-float(q01+q12-1),y_pm1+float(q01+q03));d0h=((q01+q03)<<2)+((q+q1)&3);}float x=(half_nside*(p_proj.x+p_proj.y));float y=(half_nside*(p_proj.y-p_proj.x));int i=int(x);int j=int(y);return HashDxDy((d0h<<(depth<<1))+ij2z(i,j),x-float(i),y-float(j));}uniform float opacity;vec4 get_tile_color(vec3 pos){HashDxDy result=hash_with_dxdy(0,pos.zxy);int idx=result.idx;vec2 uv=vec2(result.dy,result.dx);Tile tile=textures_tiles[idx];int idx_texture=tile.texture_idx>>6;int off=tile.texture_idx&0x3F;float idx_row=float(off>>3);float idx_col=float(off&0x7);vec2 offset=(vec2(idx_col,idx_row)+uv)*0.125;vec3 UV=vec3(offset,float(idx_texture));vec4 color=get_color_from_texture(UV);return color;}uniform sampler2D position_tex;uniform mat4 model;void main(){vec4 c=get_tile_color(normalize(frag_pos));out_frag_color=vec4(c.rgb,opacity*c.a);}`, IE = `#version 300 es
precision highp float;precision highp sampler2D;precision highp usampler2D;precision highp isampler2D;precision mediump int;in vec3 frag_pos;in vec2 out_clip_pos;out vec4 out_frag_color;struct Tile{int uniq;int texture_idx;float start_time;float empty;};uniform Tile textures_tiles[12];uniform float opacity;struct TileColor{Tile tile;vec4 color;bool found;};uniform sampler2D tex1;uniform sampler2D tex2;uniform sampler2D tex3;uniform int num_tex;uniform float scale;uniform float offset;uniform float blank;uniform float min_value;uniform float max_value;uniform int H;uniform float size_tile_uv;uniform int tex_storing_fits;uniform sampler2D colormaps;uniform float num_colormaps;uniform float colormap_id;vec4 colormap_f(float x){float id=(colormap_id+0.5)/num_colormaps;return texture(colormaps,vec2(x,id));}float linear_f(float x,float min_value,float max_value){return clamp((x-min_value)/(max_value-min_value),0.0,1.0);}float sqrt_f(float x,float min_value,float max_value){float a=linear_f(x,min_value,max_value);return sqrt(a);}float log_f(float x,float min_value,float max_value){float y=linear_f(x,min_value,max_value);float a=1000.0;return log(a*y+1.0)/log(a);}float asinh_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return asinh(10.0*d)/3.0;}float pow2_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return d*d;}float transfer_func(int H,float x,float min_value,float max_value){if(H==0){return linear_f(x,min_value,max_value);}else if(H==1){return sqrt_f(x,min_value,max_value);}else if(H==2){return log_f(x,min_value,max_value);}else if(H==3){return asinh_f(x,min_value,max_value);}else{return pow2_f(x,min_value,max_value);}}uniform float k_gamma;uniform float k_saturation;uniform float k_contrast;uniform float k_brightness;uniform float k_exposure;vec4 apply_gamma(vec4 ic,float g){float new_r=pow(ic.r,g);float new_g=pow(ic.g,g);float new_b=pow(ic.b,g);return vec4(new_r,new_g,new_b,ic.a);}vec4 apply_saturation(vec4 color,float value){const vec3 luminosity_factor=vec3(0.2126,0.7152,0.0722);vec3 grayscale=vec3(dot(color.rgb,luminosity_factor));return vec4(mix(grayscale,color.rgb,1.0+value),color.a);}vec4 apply_contrast(vec4 color,float value){return vec4(0.5+(1.0+value)*(color.rgb-0.5),color.a);}vec4 apply_brightness(vec4 color,float value){return vec4(color.rgb+value,color.a);}vec4 apply_exposure(vec4 color,float value){return vec4((1.0+value)*color.rgb,color.a);}vec4 apply_tonal(vec4 color){return apply_gamma(apply_saturation(apply_contrast(apply_brightness(color,k_brightness),k_contrast),k_saturation),k_gamma);}vec3 rgb2hsv(vec3 c){vec4 K=vec4(0.0,-1.0/3.0,2.0/3.0,-1.0);vec4 p=mix(vec4(c.bg,K.wz),vec4(c.gb,K.xy),step(c.b,c.g));vec4 q=mix(vec4(p.xyw,c.r),vec4(c.r,p.yzx),step(p.x,c.r));float d=q.x-min(q.w,q.y);float e=1.0e-10;return vec3(abs(q.z+(q.w-q.y)/(6.0*d+e)),d/(q.x+e),q.x);}vec3 hsv2rgb(vec3 c){vec4 K=vec4(1.0,2.0/3.0,1.0/3.0,3.0);vec3 p=abs(fract(c.xxx+K.xyz)*6.0-K.www);return c.z*mix(K.xxx,clamp(p-K.xxx,0.0,1.0),c.y);}vec4 get_pixels(vec3 uv){int idx_texture=int(uv.z);if(idx_texture==0){return texture(tex1,uv.xy);}else if(idx_texture==1){return texture(tex2,uv.xy);}else if(idx_texture==2){return texture(tex3,uv.xy);}else{return vec4(0.0,1.0,0.0,1.0);}}vec3 reverse_uv(vec3 uv){uv.y=size_tile_uv+2.0*size_tile_uv*floor(uv.y/size_tile_uv)-uv.y;return uv;}uniform float reversed;vec4 get_color_from_texture(vec3 UV){vec4 color=get_pixels(UV);color.r=transfer_func(H,color.r,min_value,max_value);color.g=transfer_func(H,color.g,min_value,max_value);color.b=transfer_func(H,color.b,min_value,max_value);color.rgb=mix(color.rgb,1.0-color.rgb,reversed);return apply_tonal(color);}vec4 apply_colormap_to_grayscale(float x,float a){float alpha=x*scale+offset;alpha=transfer_func(H,alpha,min_value,max_value);alpha=mix(alpha,1.0-alpha,reversed);vec4 new_color=mix(colormap_f(alpha)*a,vec4(0.0),float(x==blank||isnan(x)));return apply_tonal(new_color);}vec4 get_colormap_from_grayscale_texture(vec3 UV){vec3 uv=mix(UV,reverse_uv(UV),float(tex_storing_fits==1));vec4 color=get_pixels(uv);return apply_colormap_to_grayscale(color.r,color.a);}const float TWICE_PI=6.28318530718;const float PI=3.141592653589793;const float FOUR_OVER_PI=1.27323954474;const float TRANSITION_Z=0.66666666666;const float TRANSITION_Z_INV=1.5;int quarter(vec2 p){int x_neg=int(p.x<0.0);int y_neg=int(p.y<0.0);int q=(x_neg+y_neg)|(y_neg<<1);return q;}float xpm1(vec2 p){bool x_neg=(p.x<0.0);bool y_neg=(p.y<0.0);float lon=atan(abs(p.y),abs(p.x));float x02=lon*FOUR_OVER_PI;if(x_neg!=y_neg){return 1.0-x02;}else{return x02-1.0;}}float one_minus_z_pos(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return 1.0f-p.z;}float one_minus_z_neg(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1f){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return p.z+1.0;}int ij2z(int i,int j){int i4=i|(j<<2);int j4=(i4 ^(i4>>1))&0x22222222;int i5=i4 ^ j4 ^(j4<<1);return i5;}struct HashDxDy{int idx;float dx;float dy;};uniform sampler2D ang2pixd;HashDxDy hash_with_dxdy2(vec2 radec){vec2 aa=vec2(radec.x/TWICE_PI+1.0,(radec.y/PI)+0.5);vec3 v=texture(ang2pixd,aa).rgb;return HashDxDy(int(v.x*255.0),v.y,v.z);}HashDxDy hash_with_dxdy(int depth,vec3 p){int nside=1<<depth;float half_nside=float(nside)*0.5;float x_pm1=xpm1(p.xy);int q=quarter(p.xy);int d0h=0;vec2 p_proj=vec2(0.0);if(p.z>TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_pos(p));p_proj=vec2(x_pm1*sqrt_3_one_min_z,2.0-sqrt_3_one_min_z);d0h=q;}else if(p.z<-TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_neg(p));p_proj=vec2(x_pm1*sqrt_3_one_min_z,sqrt_3_one_min_z);d0h=q+8;}else{float y_pm1=p.z*TRANSITION_Z_INV;int q01=int(x_pm1>y_pm1);int q12=int(x_pm1>=-y_pm1);int q03=1-q12;int q1=q01&q12;p_proj=vec2(x_pm1-float(q01+q12-1),y_pm1+float(q01+q03));d0h=((q01+q03)<<2)+((q+q1)&3);}float x=(half_nside*(p_proj.x+p_proj.y));float y=(half_nside*(p_proj.y-p_proj.x));int i=int(x);int j=int(y);return HashDxDy((d0h<<(depth<<1))+ij2z(i,j),x-float(i),y-float(j));}vec4 get_tile_color(vec3 pos){HashDxDy result=hash_with_dxdy(0,pos.zxy);int idx=result.idx;vec2 uv=vec2(result.dy,result.dx);Tile tile=textures_tiles[idx];int idx_texture=tile.texture_idx>>6;int off=tile.texture_idx&0x3F;float idx_row=float(off>>3);float idx_col=float(off&0x7);vec2 offset=(vec2(idx_col,idx_row)+uv)*0.125;vec3 UV=vec3(offset,float(idx_texture));vec4 color=get_colormap_from_grayscale_texture(UV);color.a*=(1.0-tile.empty);return color;}uniform sampler2D position_tex;uniform mat4 model;void main(){vec4 c=get_tile_color(normalize(frag_pos));out_frag_color=c;out_frag_color.a=out_frag_color.a*opacity;}`, BE = `#version 300 es
precision highp float;precision highp sampler2D;precision highp usampler2D;precision highp isampler2D;precision mediump int;in vec3 frag_pos;in vec2 out_clip_pos;out vec4 out_frag_color;struct Tile{int uniq;int texture_idx;float start_time;float empty;};uniform Tile textures_tiles[12];uniform float opacity;uniform isampler2D tex1;uniform isampler2D tex2;uniform isampler2D tex3;uniform int num_tex;uniform float scale;uniform float offset;uniform float blank;uniform float min_value;uniform float max_value;uniform int H;uniform float reversed;uniform float size_tile_uv;uniform int tex_storing_fits;uniform sampler2D colormaps;uniform float num_colormaps;uniform float colormap_id;vec4 colormap_f(float x){float id=(colormap_id+0.5)/num_colormaps;return texture(colormaps,vec2(x,id));}float linear_f(float x,float min_value,float max_value){return clamp((x-min_value)/(max_value-min_value),0.0,1.0);}float sqrt_f(float x,float min_value,float max_value){float a=linear_f(x,min_value,max_value);return sqrt(a);}float log_f(float x,float min_value,float max_value){float y=linear_f(x,min_value,max_value);float a=1000.0;return log(a*y+1.0)/log(a);}float asinh_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return asinh(10.0*d)/3.0;}float pow2_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return d*d;}float transfer_func(int H,float x,float min_value,float max_value){if(H==0){return linear_f(x,min_value,max_value);}else if(H==1){return sqrt_f(x,min_value,max_value);}else if(H==2){return log_f(x,min_value,max_value);}else if(H==3){return asinh_f(x,min_value,max_value);}else{return pow2_f(x,min_value,max_value);}}uniform float k_gamma;uniform float k_saturation;uniform float k_contrast;uniform float k_brightness;uniform float k_exposure;vec4 apply_gamma(vec4 ic,float g){float new_r=pow(ic.r,g);float new_g=pow(ic.g,g);float new_b=pow(ic.b,g);return vec4(new_r,new_g,new_b,ic.a);}vec4 apply_saturation(vec4 color,float value){const vec3 luminosity_factor=vec3(0.2126,0.7152,0.0722);vec3 grayscale=vec3(dot(color.rgb,luminosity_factor));return vec4(mix(grayscale,color.rgb,1.0+value),color.a);}vec4 apply_contrast(vec4 color,float value){return vec4(0.5+(1.0+value)*(color.rgb-0.5),color.a);}vec4 apply_brightness(vec4 color,float value){return vec4(color.rgb+value,color.a);}vec4 apply_exposure(vec4 color,float value){return vec4((1.0+value)*color.rgb,color.a);}vec4 apply_tonal(vec4 color){return apply_gamma(apply_saturation(apply_contrast(apply_brightness(color,k_brightness),k_contrast),k_saturation),k_gamma);}ivec4 get_pixels(vec3 uv){int idx_texture=int(uv.z);if(idx_texture==0){return texture(tex1,uv.xy);}else if(idx_texture==1){return texture(tex2,uv.xy);}else if(idx_texture==2){return texture(tex3,uv.xy);}else{return ivec4(0,0,0,1);}}vec3 reverse_uv(vec3 uv){uv.y=size_tile_uv+2.0*size_tile_uv*floor(uv.y/size_tile_uv)-uv.y;return uv;}vec4 get_colormap_from_grayscale_texture(vec3 UV){vec3 uv=mix(UV,reverse_uv(UV),float(tex_storing_fits==1));float x=float(get_pixels(uv).r);float alpha=x*scale+offset;alpha=transfer_func(H,alpha,min_value,max_value);alpha=mix(alpha,1.0-alpha,reversed);vec4 new_color=mix(colormap_f(alpha),vec4(0.0),float(x==blank));return apply_tonal(new_color);}const float TWICE_PI=6.28318530718;const float PI=3.141592653589793;const float FOUR_OVER_PI=1.27323954474;const float TRANSITION_Z=0.66666666666;const float TRANSITION_Z_INV=1.5;int quarter(vec2 p){int x_neg=int(p.x<0.0);int y_neg=int(p.y<0.0);int q=(x_neg+y_neg)|(y_neg<<1);return q;}float xpm1(vec2 p){bool x_neg=(p.x<0.0);bool y_neg=(p.y<0.0);float lon=atan(abs(p.y),abs(p.x));float x02=lon*FOUR_OVER_PI;if(x_neg!=y_neg){return 1.0-x02;}else{return x02-1.0;}}float one_minus_z_pos(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return 1.0f-p.z;}float one_minus_z_neg(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1f){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return p.z+1.0;}int ij2z(int i,int j){int i4=i|(j<<2);int j4=(i4 ^(i4>>1))&0x22222222;int i5=i4 ^ j4 ^(j4<<1);return i5;}struct HashDxDy{int idx;float dx;float dy;};uniform sampler2D ang2pixd;HashDxDy hash_with_dxdy2(vec2 radec){vec2 aa=vec2(radec.x/TWICE_PI+1.0,(radec.y/PI)+0.5);vec3 v=texture(ang2pixd,aa).rgb;return HashDxDy(int(v.x*255.0),v.y,v.z);}HashDxDy hash_with_dxdy(int depth,vec3 p){int nside=1<<depth;float half_nside=float(nside)*0.5;float x_pm1=xpm1(p.xy);int q=quarter(p.xy);int d0h=0;vec2 p_proj=vec2(0.0);if(p.z>TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_pos(p));p_proj=vec2(x_pm1*sqrt_3_one_min_z,2.0-sqrt_3_one_min_z);d0h=q;}else if(p.z<-TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_neg(p));p_proj=vec2(x_pm1*sqrt_3_one_min_z,sqrt_3_one_min_z);d0h=q+8;}else{float y_pm1=p.z*TRANSITION_Z_INV;int q01=int(x_pm1>y_pm1);int q12=int(x_pm1>=-y_pm1);int q03=1-q12;int q1=q01&q12;p_proj=vec2(x_pm1-float(q01+q12-1),y_pm1+float(q01+q03));d0h=((q01+q03)<<2)+((q+q1)&3);}float x=(half_nside*(p_proj.x+p_proj.y));float y=(half_nside*(p_proj.y-p_proj.x));int i=int(x);int j=int(y);return HashDxDy((d0h<<(depth<<1))+ij2z(i,j),x-float(i),y-float(j));}vec4 get_tile_color(vec3 pos){HashDxDy result=hash_with_dxdy(0,pos.zxy);int idx=result.idx;vec2 uv=vec2(result.dy,result.dx);Tile tile=textures_tiles[idx];int idx_texture=tile.texture_idx>>6;int off=tile.texture_idx&0x3F;float idx_row=float(off>>3);float idx_col=float(off&0x7);vec2 offset=(vec2(idx_col,idx_row)+uv)*0.125;vec3 UV=vec3(offset,float(idx_texture));vec4 color=get_colormap_from_grayscale_texture(UV);color.a*=(1.0-tile.empty);return color;}uniform sampler2D position_tex;uniform mat4 model;void main(){vec4 c=get_tile_color(normalize(frag_pos));out_frag_color=c;out_frag_color.a=out_frag_color.a*opacity;}`, QE = `#version 300 es
precision highp float;precision highp sampler2D;precision highp usampler2D;precision highp isampler2D;precision mediump int;in vec3 frag_pos;in vec2 out_clip_pos;out vec4 out_frag_color;struct Tile{int uniq;int texture_idx;float start_time;float empty;};uniform Tile textures_tiles[12];uniform float opacity;uniform usampler2D tex1;uniform usampler2D tex2;uniform usampler2D tex3;uniform int num_tex;uniform float scale;uniform float offset;uniform float blank;uniform float min_value;uniform float max_value;uniform int H;uniform float reversed;uniform float size_tile_uv;uniform int tex_storing_fits;uniform sampler2D colormaps;uniform float num_colormaps;uniform float colormap_id;vec4 colormap_f(float x){float id=(colormap_id+0.5)/num_colormaps;return texture(colormaps,vec2(x,id));}float linear_f(float x,float min_value,float max_value){return clamp((x-min_value)/(max_value-min_value),0.0,1.0);}float sqrt_f(float x,float min_value,float max_value){float a=linear_f(x,min_value,max_value);return sqrt(a);}float log_f(float x,float min_value,float max_value){float y=linear_f(x,min_value,max_value);float a=1000.0;return log(a*y+1.0)/log(a);}float asinh_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return asinh(10.0*d)/3.0;}float pow2_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return d*d;}float transfer_func(int H,float x,float min_value,float max_value){if(H==0){return linear_f(x,min_value,max_value);}else if(H==1){return sqrt_f(x,min_value,max_value);}else if(H==2){return log_f(x,min_value,max_value);}else if(H==3){return asinh_f(x,min_value,max_value);}else{return pow2_f(x,min_value,max_value);}}uniform float k_gamma;uniform float k_saturation;uniform float k_contrast;uniform float k_brightness;uniform float k_exposure;vec4 apply_gamma(vec4 ic,float g){float new_r=pow(ic.r,g);float new_g=pow(ic.g,g);float new_b=pow(ic.b,g);return vec4(new_r,new_g,new_b,ic.a);}vec4 apply_saturation(vec4 color,float value){const vec3 luminosity_factor=vec3(0.2126,0.7152,0.0722);vec3 grayscale=vec3(dot(color.rgb,luminosity_factor));return vec4(mix(grayscale,color.rgb,1.0+value),color.a);}vec4 apply_contrast(vec4 color,float value){return vec4(0.5+(1.0+value)*(color.rgb-0.5),color.a);}vec4 apply_brightness(vec4 color,float value){return vec4(color.rgb+value,color.a);}vec4 apply_exposure(vec4 color,float value){return vec4((1.0+value)*color.rgb,color.a);}vec4 apply_tonal(vec4 color){return apply_gamma(apply_saturation(apply_contrast(apply_brightness(color,k_brightness),k_contrast),k_saturation),k_gamma);}uvec4 get_pixels(vec3 uv){int idx_texture=int(uv.z);if(idx_texture==0){return texture(tex1,uv.xy);}else if(idx_texture==1){return texture(tex2,uv.xy);}else if(idx_texture==2){return texture(tex3,uv.xy);}else{return uvec4(0,0,0,1);}}vec3 reverse_uv(vec3 uv){uv.y=size_tile_uv+2.0*size_tile_uv*floor(uv.y/size_tile_uv)-uv.y;return uv;}vec4 get_colormap_from_grayscale_texture(vec3 UV){vec3 uv=mix(UV,reverse_uv(UV),float(tex_storing_fits==1));float x=float(get_pixels(uv).r);float alpha=x*scale+offset;alpha=transfer_func(H,alpha,min_value,max_value);alpha=mix(alpha,1.0-alpha,reversed);vec4 new_color=mix(colormap_f(alpha),vec4(0.0),float(x==blank));return apply_tonal(new_color);}const float TWICE_PI=6.28318530718;const float PI=3.141592653589793;const float FOUR_OVER_PI=1.27323954474;const float TRANSITION_Z=0.66666666666;const float TRANSITION_Z_INV=1.5;int quarter(vec2 p){int x_neg=int(p.x<0.0);int y_neg=int(p.y<0.0);int q=(x_neg+y_neg)|(y_neg<<1);return q;}float xpm1(vec2 p){bool x_neg=(p.x<0.0);bool y_neg=(p.y<0.0);float lon=atan(abs(p.y),abs(p.x));float x02=lon*FOUR_OVER_PI;if(x_neg!=y_neg){return 1.0-x02;}else{return x02-1.0;}}float one_minus_z_pos(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return 1.0f-p.z;}float one_minus_z_neg(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1f){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return p.z+1.0;}int ij2z(int i,int j){int i4=i|(j<<2);int j4=(i4 ^(i4>>1))&0x22222222;int i5=i4 ^ j4 ^(j4<<1);return i5;}struct HashDxDy{int idx;float dx;float dy;};uniform sampler2D ang2pixd;HashDxDy hash_with_dxdy2(vec2 radec){vec2 aa=vec2(radec.x/TWICE_PI+1.0,(radec.y/PI)+0.5);vec3 v=texture(ang2pixd,aa).rgb;return HashDxDy(int(v.x*255.0),v.y,v.z);}HashDxDy hash_with_dxdy(int depth,vec3 p){int nside=1<<depth;float half_nside=float(nside)*0.5;float x_pm1=xpm1(p.xy);int q=quarter(p.xy);int d0h=0;vec2 p_proj=vec2(0.0);if(p.z>TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_pos(p));p_proj=vec2(x_pm1*sqrt_3_one_min_z,2.0-sqrt_3_one_min_z);d0h=q;}else if(p.z<-TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_neg(p));p_proj=vec2(x_pm1*sqrt_3_one_min_z,sqrt_3_one_min_z);d0h=q+8;}else{float y_pm1=p.z*TRANSITION_Z_INV;int q01=int(x_pm1>y_pm1);int q12=int(x_pm1>=-y_pm1);int q03=1-q12;int q1=q01&q12;p_proj=vec2(x_pm1-float(q01+q12-1),y_pm1+float(q01+q03));d0h=((q01+q03)<<2)+((q+q1)&3);}float x=(half_nside*(p_proj.x+p_proj.y));float y=(half_nside*(p_proj.y-p_proj.x));int i=int(x);int j=int(y);return HashDxDy((d0h<<(depth<<1))+ij2z(i,j),x-float(i),y-float(j));}vec4 get_tile_color(vec3 pos){HashDxDy result=hash_with_dxdy(0,pos.zxy);int idx=result.idx;vec2 uv=vec2(result.dy,result.dx);Tile tile=textures_tiles[idx];int idx_texture=tile.texture_idx>>6;int off=tile.texture_idx&0x3F;float idx_row=float(off>>3);float idx_col=float(off&0x7);vec2 offset=(vec2(idx_col,idx_row)+uv)*0.125;vec3 UV=vec3(offset,float(idx_texture));vec4 color=get_colormap_from_grayscale_texture(UV);color.a*=(1.0-tile.empty);return color;}uniform sampler2D position_tex;uniform mat4 model;void main(){vec4 c=get_tile_color(normalize(frag_pos));out_frag_color=c;out_frag_color.a=out_frag_color.a*opacity;}`, CE = `#version 300 es
precision highp float;precision mediump int;layout(location=0)in vec2 pos_clip_space;uniform vec2 ndc_to_clip;uniform float czf;void main(){gl_Position=vec4(pos_clip_space/(ndc_to_clip*czf),0.0,1.0);}`, EE = `#version 300 es
precision highp float;precision highp sampler2D;precision highp usampler2D;precision highp isampler2D;precision mediump int;out vec4 out_frag_color;uniform vec3 color;void main(){out_frag_color=vec4(color,1.0);}`, iE = `#version 300 es
precision highp float;precision mediump int;layout(location=0)in vec2 ndc_pos;layout(location=1)in vec3 uv_start;layout(location=2)in vec3 uv_end;layout(location=3)in float time_tile_received;layout(location=4)in float m0;layout(location=5)in float m1;out vec3 frag_uv_start;out vec3 frag_uv_end;out float frag_blending_factor;out float m_start;out float m_end;uniform mat4 inv_model;uniform vec2 ndc_to_clip;uniform float czf;uniform float current_time;const float PI=3.1415926535897932384626433832795;const mat4 GAL2J2000=mat4(-0.4448296299195045,0.7469822444763707,0.4941094279435681,0.0,-0.1980763734646737,0.4559837762325372,-0.8676661489811610,0.0,-0.873437090247923,-0.4838350155267381,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);const mat4 J20002GAL=mat4(-0.4448296299195045,-0.1980763734646737,-0.873437090247923,0.0,0.7469822444763707,0.4559837762325372,-0.4838350155267381,0.0,0.4941094279435681,-0.8676661489811610,-0.0548755604024359,0.0,0.0,0.0,0.0,1.0);vec2 world2clip_orthographic(vec3 p){return vec2(-p.x,p.y);}vec2 world2clip_aitoff(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float theta_by_two=theta*0.5;float alpha=acos(cos(delta)*cos(theta_by_two));float inv_sinc_alpha=1.0;if(alpha>1e-4){inv_sinc_alpha=alpha/sin(alpha);}float x=2.0*inv_sinc_alpha*cos(delta)*sin(theta_by_two);float y=inv_sinc_alpha*sin(delta);return vec2(x/PI,y/PI);}vec2 world2clip_mollweide(vec3 p){int max_iter=10;float delta=asin(p.y);float theta=atan(p.x,p.z);float cst=PI*sin(delta);float phi=delta;float dx=phi+sin(phi)-cst;int k=0;while(abs(dx)>1e-6&&k<max_iter){phi=phi-dx/(1.0+cos(phi));dx=phi+sin(phi)-cst;k=k+1;}phi=phi*0.5;float x=(-theta/PI)*cos(phi);float y=0.5*sin(phi);return vec2(x,y);}vec2 world2clip_mercator(vec3 p){float delta=asin(p.y);float theta=atan(-p.x,p.z);float x=theta/PI;float y=asinh(tan(delta/PI));return vec2(x,y);}float arc_sinc(float x){if(x>1e-4){return asin(x)/x;}else{float x2=x*x;return 1.0+x2*(1.0+x2*9.0/20.0)/6.0;}}vec2 world2clip_arc(vec3 p){if(p.z>-1.0){float r=length(p.xy);if(p.z>0.0){r=arc_sinc(r);}else{r=acos(p.z)/r;}float x=p.x*r;float y=p.y*r;return vec2(-x/PI,y/PI);}else{return vec2(1.0,0.0);}}vec2 world2clip_gnomonic(vec3 p){if(p.z<=1e-2){return vec2(1.0,0.0);}else{return vec2((-p.x/p.z)/PI,(p.y/p.z)/PI);}}const float TWICE_PI=6.28318530718;const float FOUR_OVER_PI=1.27323954474;const float TRANSITION_Z=0.66666666666;const float TRANSITION_Z_INV=1.5;float one_minus_z_pos(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return 1.0-p.z;}float one_minus_z_neg(vec3 p){float d2=dot(p.xy,p.xy);if(d2<1e-1){return d2*(0.5+d2*(0.125+d2*(0.0625+d2*(0.0390625+d2*0.02734375))));}return p.z+1.0;}vec2 xpm1_and_offset(vec2 p){int x_neg=int(p.x<0.0);int y_neg=int(p.y<0.0);int offset=-(y_neg<<2)+1+((x_neg ^ y_neg)<<1);float lon=atan(abs(p.y),abs(p.x));float x02=lon*FOUR_OVER_PI;if(x_neg!=y_neg){return vec2(1.0-x02,float(offset));}else{return vec2(x02-1.0,float(offset));}}vec2 world2clip_healpix(vec3 p){vec2 x_pm1_and_offset=xpm1_and_offset(p.xy);vec2 p_proj=vec2(0.0);if(p.z>TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_pos(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,2.0-sqrt_3_one_min_z);}else if(p.z<-TRANSITION_Z){float sqrt_3_one_min_z=sqrt(3.0*one_minus_z_neg(p));p_proj=vec2((x_pm1_and_offset.x*sqrt_3_one_min_z)+x_pm1_and_offset.y,-2.0+sqrt_3_one_min_z);}else{p_proj=vec2(atan(p.y,p.x)*FOUR_OVER_PI,p.z*TRANSITION_Z_INV);}return p_proj*vec2(-0.25,0.5);}void main(){gl_Position=vec4(ndc_pos,0.0,1.0);frag_uv_start=uv_start;frag_uv_end=uv_end;frag_blending_factor=min((current_time-time_tile_received)/200.0,1.0);m_start=m0;m_end=m1;}`, oE = `#version 300 es
precision highp float;precision highp sampler2D;precision highp isampler2D;precision mediump int;in vec3 frag_uv_start;in vec3 frag_uv_end;in float frag_blending_factor;in float m_start;in float m_end;out vec4 out_frag_color;uniform float opacity;uniform sampler2D tex1;uniform sampler2D tex2;uniform sampler2D tex3;uniform int num_tex;uniform float scale;uniform float offset;uniform float blank;uniform float min_value;uniform float max_value;uniform int H;uniform float size_tile_uv;uniform int tex_storing_fits;uniform sampler2D colormaps;uniform float num_colormaps;uniform float colormap_id;vec4 colormap_f(float x){float id=(colormap_id+0.5)/num_colormaps;return texture(colormaps,vec2(x,id));}float linear_f(float x,float min_value,float max_value){return clamp((x-min_value)/(max_value-min_value),0.0,1.0);}float sqrt_f(float x,float min_value,float max_value){float a=linear_f(x,min_value,max_value);return sqrt(a);}float log_f(float x,float min_value,float max_value){float y=linear_f(x,min_value,max_value);float a=1000.0;return log(a*y+1.0)/log(a);}float asinh_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return asinh(10.0*d)/3.0;}float pow2_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return d*d;}float transfer_func(int H,float x,float min_value,float max_value){if(H==0){return linear_f(x,min_value,max_value);}else if(H==1){return sqrt_f(x,min_value,max_value);}else if(H==2){return log_f(x,min_value,max_value);}else if(H==3){return asinh_f(x,min_value,max_value);}else{return pow2_f(x,min_value,max_value);}}uniform float k_gamma;uniform float k_saturation;uniform float k_contrast;uniform float k_brightness;uniform float k_exposure;vec4 apply_gamma(vec4 ic,float g){float new_r=pow(ic.r,g);float new_g=pow(ic.g,g);float new_b=pow(ic.b,g);return vec4(new_r,new_g,new_b,ic.a);}vec4 apply_saturation(vec4 color,float value){const vec3 luminosity_factor=vec3(0.2126,0.7152,0.0722);vec3 grayscale=vec3(dot(color.rgb,luminosity_factor));return vec4(mix(grayscale,color.rgb,1.0+value),color.a);}vec4 apply_contrast(vec4 color,float value){return vec4(0.5+(1.0+value)*(color.rgb-0.5),color.a);}vec4 apply_brightness(vec4 color,float value){return vec4(color.rgb+value,color.a);}vec4 apply_exposure(vec4 color,float value){return vec4((1.0+value)*color.rgb,color.a);}vec4 apply_tonal(vec4 color){return apply_gamma(apply_saturation(apply_contrast(apply_brightness(color,k_brightness),k_contrast),k_saturation),k_gamma);}vec3 rgb2hsv(vec3 c){vec4 K=vec4(0.0,-1.0/3.0,2.0/3.0,-1.0);vec4 p=mix(vec4(c.bg,K.wz),vec4(c.gb,K.xy),step(c.b,c.g));vec4 q=mix(vec4(p.xyw,c.r),vec4(c.r,p.yzx),step(p.x,c.r));float d=q.x-min(q.w,q.y);float e=1.0e-10;return vec3(abs(q.z+(q.w-q.y)/(6.0*d+e)),d/(q.x+e),q.x);}vec3 hsv2rgb(vec3 c){vec4 K=vec4(1.0,2.0/3.0,1.0/3.0,3.0);vec3 p=abs(fract(c.xxx+K.xyz)*6.0-K.www);return c.z*mix(K.xxx,clamp(p-K.xxx,0.0,1.0),c.y);}vec4 get_pixels(vec3 uv){int idx_texture=int(uv.z);if(idx_texture==0){return texture(tex1,uv.xy);}else if(idx_texture==1){return texture(tex2,uv.xy);}else if(idx_texture==2){return texture(tex3,uv.xy);}else{return vec4(0.0,1.0,0.0,1.0);}}vec3 reverse_uv(vec3 uv){uv.y=size_tile_uv+2.0*size_tile_uv*floor(uv.y/size_tile_uv)-uv.y;return uv;}uniform float reversed;vec4 get_color_from_texture(vec3 UV){vec4 color=get_pixels(UV);color.r=transfer_func(H,color.r,min_value,max_value);color.g=transfer_func(H,color.g,min_value,max_value);color.b=transfer_func(H,color.b,min_value,max_value);color.rgb=mix(color.rgb,1.0-color.rgb,reversed);return apply_tonal(color);}vec4 apply_colormap_to_grayscale(float x,float a){float alpha=x*scale+offset;alpha=transfer_func(H,alpha,min_value,max_value);alpha=mix(alpha,1.0-alpha,reversed);vec4 new_color=mix(colormap_f(alpha)*a,vec4(0.0),float(x==blank||isnan(x)));return apply_tonal(new_color);}vec4 get_colormap_from_grayscale_texture(vec3 UV){vec3 uv=mix(UV,reverse_uv(UV),float(tex_storing_fits==1));vec4 color=get_pixels(uv);return apply_colormap_to_grayscale(color.r,color.a);}void main(){vec4 color_start=get_color_from_texture(frag_uv_start);vec4 color_end=get_color_from_texture(frag_uv_end);out_frag_color=mix(color_start,color_end,frag_blending_factor);out_frag_color.a=opacity*out_frag_color.a;}`, DE = `#version 300 es
precision highp float;precision highp sampler2D;precision highp isampler2D;precision mediump int;in vec3 frag_uv_start;in vec3 frag_uv_end;in float frag_blending_factor;in float m_start;in float m_end;out vec4 out_frag_color;uniform sampler2D tex1;uniform sampler2D tex2;uniform sampler2D tex3;uniform int num_tex;uniform float scale;uniform float offset;uniform float blank;uniform float min_value;uniform float max_value;uniform int H;uniform float size_tile_uv;uniform int tex_storing_fits;uniform sampler2D colormaps;uniform float num_colormaps;uniform float colormap_id;vec4 colormap_f(float x){float id=(colormap_id+0.5)/num_colormaps;return texture(colormaps,vec2(x,id));}float linear_f(float x,float min_value,float max_value){return clamp((x-min_value)/(max_value-min_value),0.0,1.0);}float sqrt_f(float x,float min_value,float max_value){float a=linear_f(x,min_value,max_value);return sqrt(a);}float log_f(float x,float min_value,float max_value){float y=linear_f(x,min_value,max_value);float a=1000.0;return log(a*y+1.0)/log(a);}float asinh_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return asinh(10.0*d)/3.0;}float pow2_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return d*d;}float transfer_func(int H,float x,float min_value,float max_value){if(H==0){return linear_f(x,min_value,max_value);}else if(H==1){return sqrt_f(x,min_value,max_value);}else if(H==2){return log_f(x,min_value,max_value);}else if(H==3){return asinh_f(x,min_value,max_value);}else{return pow2_f(x,min_value,max_value);}}uniform float k_gamma;uniform float k_saturation;uniform float k_contrast;uniform float k_brightness;uniform float k_exposure;vec4 apply_gamma(vec4 ic,float g){float new_r=pow(ic.r,g);float new_g=pow(ic.g,g);float new_b=pow(ic.b,g);return vec4(new_r,new_g,new_b,ic.a);}vec4 apply_saturation(vec4 color,float value){const vec3 luminosity_factor=vec3(0.2126,0.7152,0.0722);vec3 grayscale=vec3(dot(color.rgb,luminosity_factor));return vec4(mix(grayscale,color.rgb,1.0+value),color.a);}vec4 apply_contrast(vec4 color,float value){return vec4(0.5+(1.0+value)*(color.rgb-0.5),color.a);}vec4 apply_brightness(vec4 color,float value){return vec4(color.rgb+value,color.a);}vec4 apply_exposure(vec4 color,float value){return vec4((1.0+value)*color.rgb,color.a);}vec4 apply_tonal(vec4 color){return apply_gamma(apply_saturation(apply_contrast(apply_brightness(color,k_brightness),k_contrast),k_saturation),k_gamma);}vec3 rgb2hsv(vec3 c){vec4 K=vec4(0.0,-1.0/3.0,2.0/3.0,-1.0);vec4 p=mix(vec4(c.bg,K.wz),vec4(c.gb,K.xy),step(c.b,c.g));vec4 q=mix(vec4(p.xyw,c.r),vec4(c.r,p.yzx),step(p.x,c.r));float d=q.x-min(q.w,q.y);float e=1.0e-10;return vec3(abs(q.z+(q.w-q.y)/(6.0*d+e)),d/(q.x+e),q.x);}vec3 hsv2rgb(vec3 c){vec4 K=vec4(1.0,2.0/3.0,1.0/3.0,3.0);vec3 p=abs(fract(c.xxx+K.xyz)*6.0-K.www);return c.z*mix(K.xxx,clamp(p-K.xxx,0.0,1.0),c.y);}vec4 get_pixels(vec3 uv){int idx_texture=int(uv.z);if(idx_texture==0){return texture(tex1,uv.xy);}else if(idx_texture==1){return texture(tex2,uv.xy);}else if(idx_texture==2){return texture(tex3,uv.xy);}else{return vec4(0.0,1.0,0.0,1.0);}}vec3 reverse_uv(vec3 uv){uv.y=size_tile_uv+2.0*size_tile_uv*floor(uv.y/size_tile_uv)-uv.y;return uv;}uniform float reversed;vec4 get_color_from_texture(vec3 UV){vec4 color=get_pixels(UV);color.r=transfer_func(H,color.r,min_value,max_value);color.g=transfer_func(H,color.g,min_value,max_value);color.b=transfer_func(H,color.b,min_value,max_value);color.rgb=mix(color.rgb,1.0-color.rgb,reversed);return apply_tonal(color);}vec4 apply_colormap_to_grayscale(float x,float a){float alpha=x*scale+offset;alpha=transfer_func(H,alpha,min_value,max_value);alpha=mix(alpha,1.0-alpha,reversed);vec4 new_color=mix(colormap_f(alpha)*a,vec4(0.0),float(x==blank||isnan(x)));return apply_tonal(new_color);}vec4 get_colormap_from_grayscale_texture(vec3 UV){vec3 uv=mix(UV,reverse_uv(UV),float(tex_storing_fits==1));vec4 color=get_pixels(uv);return apply_colormap_to_grayscale(color.r,color.a);}uniform float opacity;void main(){vec4 color_start=get_colormap_from_grayscale_texture(frag_uv_start);vec4 color_end=get_colormap_from_grayscale_texture(frag_uv_end);out_frag_color=mix(color_start,color_end,frag_blending_factor);out_frag_color.a=out_frag_color.a*opacity;}`, aE = `#version 300 es
precision highp float;precision highp sampler2D;precision highp isampler2D;precision mediump int;in vec3 frag_uv_start;in vec3 frag_uv_end;in float frag_blending_factor;in float m_start;in float m_end;out vec4 out_frag_color;uniform isampler2D tex1;uniform isampler2D tex2;uniform isampler2D tex3;uniform int num_tex;uniform float scale;uniform float offset;uniform float blank;uniform float min_value;uniform float max_value;uniform int H;uniform float reversed;uniform float size_tile_uv;uniform int tex_storing_fits;uniform sampler2D colormaps;uniform float num_colormaps;uniform float colormap_id;vec4 colormap_f(float x){float id=(colormap_id+0.5)/num_colormaps;return texture(colormaps,vec2(x,id));}float linear_f(float x,float min_value,float max_value){return clamp((x-min_value)/(max_value-min_value),0.0,1.0);}float sqrt_f(float x,float min_value,float max_value){float a=linear_f(x,min_value,max_value);return sqrt(a);}float log_f(float x,float min_value,float max_value){float y=linear_f(x,min_value,max_value);float a=1000.0;return log(a*y+1.0)/log(a);}float asinh_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return asinh(10.0*d)/3.0;}float pow2_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return d*d;}float transfer_func(int H,float x,float min_value,float max_value){if(H==0){return linear_f(x,min_value,max_value);}else if(H==1){return sqrt_f(x,min_value,max_value);}else if(H==2){return log_f(x,min_value,max_value);}else if(H==3){return asinh_f(x,min_value,max_value);}else{return pow2_f(x,min_value,max_value);}}uniform float k_gamma;uniform float k_saturation;uniform float k_contrast;uniform float k_brightness;uniform float k_exposure;vec4 apply_gamma(vec4 ic,float g){float new_r=pow(ic.r,g);float new_g=pow(ic.g,g);float new_b=pow(ic.b,g);return vec4(new_r,new_g,new_b,ic.a);}vec4 apply_saturation(vec4 color,float value){const vec3 luminosity_factor=vec3(0.2126,0.7152,0.0722);vec3 grayscale=vec3(dot(color.rgb,luminosity_factor));return vec4(mix(grayscale,color.rgb,1.0+value),color.a);}vec4 apply_contrast(vec4 color,float value){return vec4(0.5+(1.0+value)*(color.rgb-0.5),color.a);}vec4 apply_brightness(vec4 color,float value){return vec4(color.rgb+value,color.a);}vec4 apply_exposure(vec4 color,float value){return vec4((1.0+value)*color.rgb,color.a);}vec4 apply_tonal(vec4 color){return apply_gamma(apply_saturation(apply_contrast(apply_brightness(color,k_brightness),k_contrast),k_saturation),k_gamma);}ivec4 get_pixels(vec3 uv){int idx_texture=int(uv.z);if(idx_texture==0){return texture(tex1,uv.xy);}else if(idx_texture==1){return texture(tex2,uv.xy);}else if(idx_texture==2){return texture(tex3,uv.xy);}else{return ivec4(0,0,0,1);}}vec3 reverse_uv(vec3 uv){uv.y=size_tile_uv+2.0*size_tile_uv*floor(uv.y/size_tile_uv)-uv.y;return uv;}vec4 get_colormap_from_grayscale_texture(vec3 UV){vec3 uv=mix(UV,reverse_uv(UV),float(tex_storing_fits==1));float x=float(get_pixels(uv).r);float alpha=x*scale+offset;alpha=transfer_func(H,alpha,min_value,max_value);alpha=mix(alpha,1.0-alpha,reversed);vec4 new_color=mix(colormap_f(alpha),vec4(0.0),float(x==blank));return apply_tonal(new_color);}uniform float opacity;void main(){vec4 color_start=get_colormap_from_grayscale_texture(frag_uv_start);vec4 color_end=get_colormap_from_grayscale_texture(frag_uv_end);out_frag_color=mix(color_start,color_end,frag_blending_factor);out_frag_color.a=out_frag_color.a*opacity;}`, wE = `#version 300 es
precision highp float;precision highp sampler2D;precision highp isampler2D;precision highp usampler2D;precision mediump int;in vec3 frag_uv_start;in vec3 frag_uv_end;in float frag_blending_factor;in float m_start;in float m_end;out vec4 out_frag_color;uniform usampler2D tex1;uniform usampler2D tex2;uniform usampler2D tex3;uniform int num_tex;uniform float scale;uniform float offset;uniform float blank;uniform float min_value;uniform float max_value;uniform int H;uniform float reversed;uniform float size_tile_uv;uniform int tex_storing_fits;uniform sampler2D colormaps;uniform float num_colormaps;uniform float colormap_id;vec4 colormap_f(float x){float id=(colormap_id+0.5)/num_colormaps;return texture(colormaps,vec2(x,id));}float linear_f(float x,float min_value,float max_value){return clamp((x-min_value)/(max_value-min_value),0.0,1.0);}float sqrt_f(float x,float min_value,float max_value){float a=linear_f(x,min_value,max_value);return sqrt(a);}float log_f(float x,float min_value,float max_value){float y=linear_f(x,min_value,max_value);float a=1000.0;return log(a*y+1.0)/log(a);}float asinh_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return asinh(10.0*d)/3.0;}float pow2_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return d*d;}float transfer_func(int H,float x,float min_value,float max_value){if(H==0){return linear_f(x,min_value,max_value);}else if(H==1){return sqrt_f(x,min_value,max_value);}else if(H==2){return log_f(x,min_value,max_value);}else if(H==3){return asinh_f(x,min_value,max_value);}else{return pow2_f(x,min_value,max_value);}}uniform float k_gamma;uniform float k_saturation;uniform float k_contrast;uniform float k_brightness;uniform float k_exposure;vec4 apply_gamma(vec4 ic,float g){float new_r=pow(ic.r,g);float new_g=pow(ic.g,g);float new_b=pow(ic.b,g);return vec4(new_r,new_g,new_b,ic.a);}vec4 apply_saturation(vec4 color,float value){const vec3 luminosity_factor=vec3(0.2126,0.7152,0.0722);vec3 grayscale=vec3(dot(color.rgb,luminosity_factor));return vec4(mix(grayscale,color.rgb,1.0+value),color.a);}vec4 apply_contrast(vec4 color,float value){return vec4(0.5+(1.0+value)*(color.rgb-0.5),color.a);}vec4 apply_brightness(vec4 color,float value){return vec4(color.rgb+value,color.a);}vec4 apply_exposure(vec4 color,float value){return vec4((1.0+value)*color.rgb,color.a);}vec4 apply_tonal(vec4 color){return apply_gamma(apply_saturation(apply_contrast(apply_brightness(color,k_brightness),k_contrast),k_saturation),k_gamma);}uvec4 get_pixels(vec3 uv){int idx_texture=int(uv.z);if(idx_texture==0){return texture(tex1,uv.xy);}else if(idx_texture==1){return texture(tex2,uv.xy);}else if(idx_texture==2){return texture(tex3,uv.xy);}else{return uvec4(0,0,0,1);}}vec3 reverse_uv(vec3 uv){uv.y=size_tile_uv+2.0*size_tile_uv*floor(uv.y/size_tile_uv)-uv.y;return uv;}vec4 get_colormap_from_grayscale_texture(vec3 UV){vec3 uv=mix(UV,reverse_uv(UV),float(tex_storing_fits==1));float x=float(get_pixels(uv).r);float alpha=x*scale+offset;alpha=transfer_func(H,alpha,min_value,max_value);alpha=mix(alpha,1.0-alpha,reversed);vec4 new_color=mix(colormap_f(alpha),vec4(0.0),float(x==blank));return apply_tonal(new_color);}uniform float opacity;void main(){vec4 color_start=get_colormap_from_grayscale_texture(frag_uv_start);vec4 color_end=get_colormap_from_grayscale_texture(frag_uv_end);out_frag_color=mix(color_start,color_end,frag_blending_factor);out_frag_color.a=out_frag_color.a*opacity;}`, tE = `#version 300 es
precision mediump float;layout(location=0)in vec2 a_pos;out vec2 v_tc;void main(){gl_Position=vec4(a_pos*2.-1.,0.0,1.0);v_tc=a_pos;}`, sE = `#version 300 es
precision mediump float;in vec2 v_tc;out vec4 color;uniform sampler2D fbo_tex;vec3 srgb_from_linear(vec3 rgb){bvec3 cutoff=lessThan(rgb,vec3(0.0031308));vec3 lower=rgb*vec3(3294.6);vec3 higher=vec3(269.025)*pow(rgb,vec3(1.0/2.4))-vec3(14.025);return mix(higher,lower,vec3(cutoff));}vec4 srgba_from_linear(vec4 rgba){return vec4(srgb_from_linear(rgba.rgb),255.0*rgba.a);}void main(){color=texture(fbo_tex,v_tc);}`, eE = `#version 300 es
precision highp float;precision mediump int;layout(location=0)in vec2 ndc_pos;layout(location=1)in vec2 uv;out vec2 frag_uv;void main(){gl_Position=vec4(ndc_pos,0.0,1.0);frag_uv=uv;}`, cE = `#version 300 es
precision highp float;precision highp sampler2D;precision highp isampler2D;precision highp usampler2D;precision mediump int;out vec4 out_frag_color;in vec2 frag_uv;uniform sampler2D tex;uniform float opacity;uniform sampler2D tex1;uniform sampler2D tex2;uniform sampler2D tex3;uniform int num_tex;uniform float scale;uniform float offset;uniform float blank;uniform float min_value;uniform float max_value;uniform int H;uniform float size_tile_uv;uniform int tex_storing_fits;uniform sampler2D colormaps;uniform float num_colormaps;uniform float colormap_id;vec4 colormap_f(float x){float id=(colormap_id+0.5)/num_colormaps;return texture(colormaps,vec2(x,id));}float linear_f(float x,float min_value,float max_value){return clamp((x-min_value)/(max_value-min_value),0.0,1.0);}float sqrt_f(float x,float min_value,float max_value){float a=linear_f(x,min_value,max_value);return sqrt(a);}float log_f(float x,float min_value,float max_value){float y=linear_f(x,min_value,max_value);float a=1000.0;return log(a*y+1.0)/log(a);}float asinh_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return asinh(10.0*d)/3.0;}float pow2_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return d*d;}float transfer_func(int H,float x,float min_value,float max_value){if(H==0){return linear_f(x,min_value,max_value);}else if(H==1){return sqrt_f(x,min_value,max_value);}else if(H==2){return log_f(x,min_value,max_value);}else if(H==3){return asinh_f(x,min_value,max_value);}else{return pow2_f(x,min_value,max_value);}}uniform float k_gamma;uniform float k_saturation;uniform float k_contrast;uniform float k_brightness;uniform float k_exposure;vec4 apply_gamma(vec4 ic,float g){float new_r=pow(ic.r,g);float new_g=pow(ic.g,g);float new_b=pow(ic.b,g);return vec4(new_r,new_g,new_b,ic.a);}vec4 apply_saturation(vec4 color,float value){const vec3 luminosity_factor=vec3(0.2126,0.7152,0.0722);vec3 grayscale=vec3(dot(color.rgb,luminosity_factor));return vec4(mix(grayscale,color.rgb,1.0+value),color.a);}vec4 apply_contrast(vec4 color,float value){return vec4(0.5+(1.0+value)*(color.rgb-0.5),color.a);}vec4 apply_brightness(vec4 color,float value){return vec4(color.rgb+value,color.a);}vec4 apply_exposure(vec4 color,float value){return vec4((1.0+value)*color.rgb,color.a);}vec4 apply_tonal(vec4 color){return apply_gamma(apply_saturation(apply_contrast(apply_brightness(color,k_brightness),k_contrast),k_saturation),k_gamma);}vec3 rgb2hsv(vec3 c){vec4 K=vec4(0.0,-1.0/3.0,2.0/3.0,-1.0);vec4 p=mix(vec4(c.bg,K.wz),vec4(c.gb,K.xy),step(c.b,c.g));vec4 q=mix(vec4(p.xyw,c.r),vec4(c.r,p.yzx),step(p.x,c.r));float d=q.x-min(q.w,q.y);float e=1.0e-10;return vec3(abs(q.z+(q.w-q.y)/(6.0*d+e)),d/(q.x+e),q.x);}vec3 hsv2rgb(vec3 c){vec4 K=vec4(1.0,2.0/3.0,1.0/3.0,3.0);vec3 p=abs(fract(c.xxx+K.xyz)*6.0-K.www);return c.z*mix(K.xxx,clamp(p-K.xxx,0.0,1.0),c.y);}vec4 get_pixels(vec3 uv){int idx_texture=int(uv.z);if(idx_texture==0){return texture(tex1,uv.xy);}else if(idx_texture==1){return texture(tex2,uv.xy);}else if(idx_texture==2){return texture(tex3,uv.xy);}else{return vec4(0.0,1.0,0.0,1.0);}}vec3 reverse_uv(vec3 uv){uv.y=size_tile_uv+2.0*size_tile_uv*floor(uv.y/size_tile_uv)-uv.y;return uv;}uniform float reversed;vec4 get_color_from_texture(vec3 UV){vec4 color=get_pixels(UV);color.r=transfer_func(H,color.r,min_value,max_value);color.g=transfer_func(H,color.g,min_value,max_value);color.b=transfer_func(H,color.b,min_value,max_value);color.rgb=mix(color.rgb,1.0-color.rgb,reversed);return apply_tonal(color);}vec4 apply_colormap_to_grayscale(float x,float a){float alpha=x*scale+offset;alpha=transfer_func(H,alpha,min_value,max_value);alpha=mix(alpha,1.0-alpha,reversed);vec4 new_color=mix(colormap_f(alpha)*a,vec4(0.0),float(x==blank||isnan(x)));return apply_tonal(new_color);}vec4 get_colormap_from_grayscale_texture(vec3 UV){vec3 uv=mix(UV,reverse_uv(UV),float(tex_storing_fits==1));vec4 color=get_pixels(uv);return apply_colormap_to_grayscale(color.r,color.a);}void main(){vec4 color=texture(tex,frag_uv);out_frag_color=apply_colormap_to_grayscale(color.r,color.a);out_frag_color.a=out_frag_color.a*opacity;}`, GE = `#version 300 es
precision highp float;precision highp sampler2D;precision highp isampler2D;precision highp usampler2D;precision mediump int;out vec4 out_frag_color;in vec2 frag_uv;uniform usampler2D tex;uniform float opacity;uniform sampler2D tex1;uniform sampler2D tex2;uniform sampler2D tex3;uniform int num_tex;uniform float scale;uniform float offset;uniform float blank;uniform float min_value;uniform float max_value;uniform int H;uniform float size_tile_uv;uniform int tex_storing_fits;uniform sampler2D colormaps;uniform float num_colormaps;uniform float colormap_id;vec4 colormap_f(float x){float id=(colormap_id+0.5)/num_colormaps;return texture(colormaps,vec2(x,id));}float linear_f(float x,float min_value,float max_value){return clamp((x-min_value)/(max_value-min_value),0.0,1.0);}float sqrt_f(float x,float min_value,float max_value){float a=linear_f(x,min_value,max_value);return sqrt(a);}float log_f(float x,float min_value,float max_value){float y=linear_f(x,min_value,max_value);float a=1000.0;return log(a*y+1.0)/log(a);}float asinh_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return asinh(10.0*d)/3.0;}float pow2_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return d*d;}float transfer_func(int H,float x,float min_value,float max_value){if(H==0){return linear_f(x,min_value,max_value);}else if(H==1){return sqrt_f(x,min_value,max_value);}else if(H==2){return log_f(x,min_value,max_value);}else if(H==3){return asinh_f(x,min_value,max_value);}else{return pow2_f(x,min_value,max_value);}}uniform float k_gamma;uniform float k_saturation;uniform float k_contrast;uniform float k_brightness;uniform float k_exposure;vec4 apply_gamma(vec4 ic,float g){float new_r=pow(ic.r,g);float new_g=pow(ic.g,g);float new_b=pow(ic.b,g);return vec4(new_r,new_g,new_b,ic.a);}vec4 apply_saturation(vec4 color,float value){const vec3 luminosity_factor=vec3(0.2126,0.7152,0.0722);vec3 grayscale=vec3(dot(color.rgb,luminosity_factor));return vec4(mix(grayscale,color.rgb,1.0+value),color.a);}vec4 apply_contrast(vec4 color,float value){return vec4(0.5+(1.0+value)*(color.rgb-0.5),color.a);}vec4 apply_brightness(vec4 color,float value){return vec4(color.rgb+value,color.a);}vec4 apply_exposure(vec4 color,float value){return vec4((1.0+value)*color.rgb,color.a);}vec4 apply_tonal(vec4 color){return apply_gamma(apply_saturation(apply_contrast(apply_brightness(color,k_brightness),k_contrast),k_saturation),k_gamma);}vec3 rgb2hsv(vec3 c){vec4 K=vec4(0.0,-1.0/3.0,2.0/3.0,-1.0);vec4 p=mix(vec4(c.bg,K.wz),vec4(c.gb,K.xy),step(c.b,c.g));vec4 q=mix(vec4(p.xyw,c.r),vec4(c.r,p.yzx),step(p.x,c.r));float d=q.x-min(q.w,q.y);float e=1.0e-10;return vec3(abs(q.z+(q.w-q.y)/(6.0*d+e)),d/(q.x+e),q.x);}vec3 hsv2rgb(vec3 c){vec4 K=vec4(1.0,2.0/3.0,1.0/3.0,3.0);vec3 p=abs(fract(c.xxx+K.xyz)*6.0-K.www);return c.z*mix(K.xxx,clamp(p-K.xxx,0.0,1.0),c.y);}vec4 get_pixels(vec3 uv){int idx_texture=int(uv.z);if(idx_texture==0){return texture(tex1,uv.xy);}else if(idx_texture==1){return texture(tex2,uv.xy);}else if(idx_texture==2){return texture(tex3,uv.xy);}else{return vec4(0.0,1.0,0.0,1.0);}}vec3 reverse_uv(vec3 uv){uv.y=size_tile_uv+2.0*size_tile_uv*floor(uv.y/size_tile_uv)-uv.y;return uv;}uniform float reversed;vec4 get_color_from_texture(vec3 UV){vec4 color=get_pixels(UV);color.r=transfer_func(H,color.r,min_value,max_value);color.g=transfer_func(H,color.g,min_value,max_value);color.b=transfer_func(H,color.b,min_value,max_value);color.rgb=mix(color.rgb,1.0-color.rgb,reversed);return apply_tonal(color);}vec4 apply_colormap_to_grayscale(float x,float a){float alpha=x*scale+offset;alpha=transfer_func(H,alpha,min_value,max_value);alpha=mix(alpha,1.0-alpha,reversed);vec4 new_color=mix(colormap_f(alpha)*a,vec4(0.0),float(x==blank||isnan(x)));return apply_tonal(new_color);}vec4 get_colormap_from_grayscale_texture(vec3 UV){vec3 uv=mix(UV,reverse_uv(UV),float(tex_storing_fits==1));vec4 color=get_pixels(uv);return apply_colormap_to_grayscale(color.r,color.a);}void main(){uvec4 color=texture(tex,frag_uv);out_frag_color=apply_colormap_to_grayscale(float(color.r),float(color.a));out_frag_color.a=out_frag_color.a*opacity;}`, hE = `#version 300 es
precision highp float;precision highp sampler2D;precision highp isampler2D;precision highp usampler2D;precision mediump int;out vec4 out_frag_color;in vec2 frag_uv;uniform isampler2D tex;uniform float opacity;uniform sampler2D tex1;uniform sampler2D tex2;uniform sampler2D tex3;uniform int num_tex;uniform float scale;uniform float offset;uniform float blank;uniform float min_value;uniform float max_value;uniform int H;uniform float size_tile_uv;uniform int tex_storing_fits;uniform sampler2D colormaps;uniform float num_colormaps;uniform float colormap_id;vec4 colormap_f(float x){float id=(colormap_id+0.5)/num_colormaps;return texture(colormaps,vec2(x,id));}float linear_f(float x,float min_value,float max_value){return clamp((x-min_value)/(max_value-min_value),0.0,1.0);}float sqrt_f(float x,float min_value,float max_value){float a=linear_f(x,min_value,max_value);return sqrt(a);}float log_f(float x,float min_value,float max_value){float y=linear_f(x,min_value,max_value);float a=1000.0;return log(a*y+1.0)/log(a);}float asinh_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return asinh(10.0*d)/3.0;}float pow2_f(float x,float min_value,float max_value){float d=linear_f(x,min_value,max_value);return d*d;}float transfer_func(int H,float x,float min_value,float max_value){if(H==0){return linear_f(x,min_value,max_value);}else if(H==1){return sqrt_f(x,min_value,max_value);}else if(H==2){return log_f(x,min_value,max_value);}else if(H==3){return asinh_f(x,min_value,max_value);}else{return pow2_f(x,min_value,max_value);}}uniform float k_gamma;uniform float k_saturation;uniform float k_contrast;uniform float k_brightness;uniform float k_exposure;vec4 apply_gamma(vec4 ic,float g){float new_r=pow(ic.r,g);float new_g=pow(ic.g,g);float new_b=pow(ic.b,g);return vec4(new_r,new_g,new_b,ic.a);}vec4 apply_saturation(vec4 color,float value){const vec3 luminosity_factor=vec3(0.2126,0.7152,0.0722);vec3 grayscale=vec3(dot(color.rgb,luminosity_factor));return vec4(mix(grayscale,color.rgb,1.0+value),color.a);}vec4 apply_contrast(vec4 color,float value){return vec4(0.5+(1.0+value)*(color.rgb-0.5),color.a);}vec4 apply_brightness(vec4 color,float value){return vec4(color.rgb+value,color.a);}vec4 apply_exposure(vec4 color,float value){return vec4((1.0+value)*color.rgb,color.a);}vec4 apply_tonal(vec4 color){return apply_gamma(apply_saturation(apply_contrast(apply_brightness(color,k_brightness),k_contrast),k_saturation),k_gamma);}vec3 rgb2hsv(vec3 c){vec4 K=vec4(0.0,-1.0/3.0,2.0/3.0,-1.0);vec4 p=mix(vec4(c.bg,K.wz),vec4(c.gb,K.xy),step(c.b,c.g));vec4 q=mix(vec4(p.xyw,c.r),vec4(c.r,p.yzx),step(p.x,c.r));float d=q.x-min(q.w,q.y);float e=1.0e-10;return vec3(abs(q.z+(q.w-q.y)/(6.0*d+e)),d/(q.x+e),q.x);}vec3 hsv2rgb(vec3 c){vec4 K=vec4(1.0,2.0/3.0,1.0/3.0,3.0);vec3 p=abs(fract(c.xxx+K.xyz)*6.0-K.www);return c.z*mix(K.xxx,clamp(p-K.xxx,0.0,1.0),c.y);}vec4 get_pixels(vec3 uv){int idx_texture=int(uv.z);if(idx_texture==0){return texture(tex1,uv.xy);}else if(idx_texture==1){return texture(tex2,uv.xy);}else if(idx_texture==2){return texture(tex3,uv.xy);}else{return vec4(0.0,1.0,0.0,1.0);}}vec3 reverse_uv(vec3 uv){uv.y=size_tile_uv+2.0*size_tile_uv*floor(uv.y/size_tile_uv)-uv.y;return uv;}uniform float reversed;vec4 get_color_from_texture(vec3 UV){vec4 color=get_pixels(UV);color.r=transfer_func(H,color.r,min_value,max_value);color.g=transfer_func(H,color.g,min_value,max_value);color.b=transfer_func(H,color.b,min_value,max_value);color.rgb=mix(color.rgb,1.0-color.rgb,reversed);return apply_tonal(color);}vec4 apply_colormap_to_grayscale(float x,float a){float alpha=x*scale+offset;alpha=transfer_func(H,alpha,min_value,max_value);alpha=mix(alpha,1.0-alpha,reversed);vec4 new_color=mix(colormap_f(alpha)*a,vec4(0.0),float(x==blank||isnan(x)));return apply_tonal(new_color);}vec4 get_colormap_from_grayscale_texture(vec3 UV){vec3 uv=mix(UV,reverse_uv(UV),float(tex_storing_fits==1));vec4 color=get_pixels(uv);return apply_colormap_to_grayscale(color.r,color.a);}void main(){ivec4 color=texture(tex,frag_uv);out_frag_color=apply_colormap_to_grayscale(float(color.r),float(color.a));out_frag_color.a=out_frag_color.a*opacity;}`;
let ME = [
  // Catalog shaders
  {
    id: "CatalogAitoffVS",
    content: OC
  },
  {
    id: "CatalogHEALPixVS",
    content: WC
  },
  {
    id: "CatalogMercatVS",
    content: vC
  },
  {
    id: "CatalogArcVS",
    content: bC
  },
  {
    id: "CatalogTanVS",
    content: ZC
  },
  {
    id: "CatalogMollVS",
    content: TC
  },
  {
    id: "CatalogOrthoVS",
    content: VC
  },
  {
    id: "CatalogOrthoFS",
    content: jC
  },
  {
    id: "CatalogFS",
    content: _C
  },
  // Colormap shaders
  {
    id: "ColormapCatalogVS",
    content: PC
  },
  {
    id: "ColormapCatalogFS",
    content: XC
  },
  // Grid shader
  {
    id: "GridVS_CPU",
    content: zC
  },
  {
    id: "GridFS_CPU",
    content: $C
  },
  // HiPS shaders
  // Raytracer
  {
    id: "RayTracerVS",
    content: AE
  },
  {
    id: "RayTracerColorFS",
    content: gE
  },
  {
    id: "RayTracerGrayscale2ColormapFS",
    content: IE
  },
  {
    id: "RayTracerGrayscale2ColormapIntegerFS",
    content: BE
  },
  {
    id: "RayTracerGrayscale2ColormapUnsignedFS",
    content: QE
  },
  {
    id: "RayTracerFontVS",
    content: CE
  },
  {
    id: "RayTracerFontFS",
    content: EE
  },
  /// Rasterizer
  {
    id: "RasterizerVS",
    content: iE
  },
  {
    id: "RasterizerColorFS",
    content: oE
  },
  {
    id: "RasterizerGrayscale2ColormapFS",
    content: DE
  },
  {
    id: "RasterizerGrayscale2ColormapIntegerFS",
    content: aE
  },
  {
    id: "RasterizerGrayscale2ColormapUnsignedFS",
    content: wE
  },
  // Post
  {
    id: "PostVS",
    content: tE
  },
  {
    id: "PostFS",
    content: sE
  },
  // Fits
  {
    id: "FitsVS",
    content: eE
  },
  {
    id: "FitsFS",
    content: cE
  },
  {
    id: "FitsFSUnsigned",
    content: GE
  },
  {
    id: "FitsFSInteger",
    content: hE
  }
];
function yE() {
  return ME;
}
const rE = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACAAAAAgCAYAAABzenr0AAAABmJLR0QA/wD/AP+gvaeTAAAACXBIWXMAAAsTAAALEwEAmpwYAAAAB3RJTUUH4wwFDTQFbbJKlAAAAB1pVFh0Q29tbWVudAAAAAAAQ3JlYXRlZCB3aXRoIEdJTVBkLmUHAAAC+klEQVRYw8WX627bMAyFD23ZTppi2Ps/5YAt8UWKuR8jvRNWSrd2wAoQgRtH/HhIURTwn//kb15W1Q7AAKA36+jrHcDdLIvI/s8AVHUAMJoNAJIBiJma3QEUABnABmATkfxhAIv4DOAEYCKAwaLv/VUCyASwAlgAzC1F0jtRvxiAQ7gKySwqUMw2s8XfV9VbTY3UcD4CuBjAJQBMBNARwE4AqwGMlK5OVa8isj0FsMgvAF4J4myfkxkXopD8noLVLIX3oKrKSqRKzj3qaLVa6AKAF+Bi1lOqfKfsqvrdayIqcCa7NJQ4EUBUIJP0A21Vrah0ewAw6U8NCAfx/7sKDFBI/oFq5IicaiSr6iYihRUYKcdngnkhFVyJExWYUBPyyk8UOTco3x3u5xeA5X4KECezF6qL15CKgQAKRS/UH7g3eHGOACZVnRMVYwrdbgowEcJV6AwgA5hDc3KoiRwPpN4QAdyGBgSnxLdlbwALOeftGDvoQLujd4CeLFVspPScQj0ki9TroVC07DgFHz0DeEvtKuZbiWEmc/7FnldbJ1cOrL6y1uGvCwAIzwwWAR3mK7VndiYVi2s/nOcPTbFhe9jTK4Bv9lks9/59aw32cTQiXvheeb6Hk24xmWGVX+zT97k3pUK/v4f1FcCeqGqjlcbZPhJ4pm24WHtd6DTMDZDDGKBwq6w4X6nL+W+WSh+4ma0EstFaD+pwJA+jVIg49nZvqwOdBb7vb2QzKbLS+keACQBEZFfVjYifOd9pu6VwGG2UigixxPSIyM6H0frEcZx61ndOQ0/FlUBYCffz+zgWkayqS6OJxDN/ejKQLA2ImRURkVIbSGbqeH2YZGJ/jyPZXpuGCeKHPbu9HcmsFm6VBsXnOU/GtfTkAMA74wrgxiP6m6HUUnENzgtJH8ctCYNHpmJkiCuAaxzNq2O5iGyqGkftGH0LoITeccgeR/KP3Iwmmhfi1WynOuBa8KLbP3M3TEH+VLmc3kMdfP5u2FAkPbmcHgB/ejv+CcL32flvuhUBAAAAAElFTkSuQmCC";
let NE = function() {
  function B(A, g) {
    const C = yE();
    this.webclient = new A.WebClient(
      g,
      C,
      {
        kernel: rE
      }
    );
  }
  return B;
}(), cB = {};
cB.log = function(B, A) {
  try {
    var g = "//alasky.unistra.fr/cgi/AladinLiteLogger/log.py", C = "";
    A && (C = JSON.stringify(A)), O.ajax({
      url: g,
      data: { action: B, params: C, pageUrl: window.location.href, referer: document.referrer ? document.referrer : "" },
      method: "GET",
      dataType: "json"
      // as alasky supports CORS, we do not need JSONP any longer
    });
  } catch (I) {
    window.console && console.log("Exception: " + I);
  }
};
class CA {
  static AL_USE_WASM = new CA("AL:Wasm");
  static LOADING_STATE = new CA("AL:Layer.loading");
  static BACKGROUND_COLOR_CHANGED = new CA("AL:BackgroundColor.changed");
  static COO_GRID_ENABLED = new CA("AL:cooGrid.enabled");
  static COO_GRID_DISABLED = new CA("AL:cooGrid.disabled");
  static COO_GRID_UPDATED = new CA("AL:cooGrid.updated");
  static PROJECTION_CHANGED = new CA("AL:projection.changed");
  static HIPS_LAYER_ADDED = new CA("AL:HiPSLayer.added");
  static HIPS_LAYER_REMOVED = new CA("AL:HiPSLayer.removed");
  static HIPS_LAYER_RENAMED = new CA("AL:HiPSLayer.renamed");
  static HIPS_LAYER_SWAP = new CA("AL:HiPSLayer.swap");
  static HIPS_LAYER_CHANGED = new CA("AL:HiPSLayer.changed");
  static GRAPHIC_OVERLAY_LAYER_ADDED = new CA("AL:GraphicOverlayLayer.added");
  static GRAPHIC_OVERLAY_LAYER_REMOVED = new CA("AL:GraphicOverlayLayer.removed");
  static GRAPHIC_OVERLAY_LAYER_CHANGED = new CA("AL:GraphicOverlayLayer.changed");
  constructor(A) {
    this.name = A;
  }
  dispatchedTo(A, g) {
    g ? A.dispatchEvent(new CustomEvent(this.name, { detail: g })) : A.dispatchEvent(new CustomEvent(this.name));
  }
  listenedBy(A, g) {
    A.addEventListener(this.name, g);
  }
  remove(A, g) {
    A.removeEventListener(this.name, g);
  }
}
let OI = function() {
  function B(g) {
    this.opacity = g && g.opacity || 1, this.colormap = g && g.colormap || "native", this.colormap = A(this.colormap), this.stretch = g && g.stretch || "linear", this.stretch = this.stretch.toLowerCase(), this.reversed = !1, g && g.reversed === !0 && (this.reversed = !0), this.minCut = g && g.minCut || 0, this.maxCut = g && g.maxCut || 1, this.additiveBlending = g && g.additive, this.additiveBlending === void 0 && (this.additiveBlending = !1), this.kGamma = g && g.gamma || 1, this.kSaturation = g && g.saturation || 0, this.kBrightness = g && g.brightness || 0, this.kContrast = g && g.contrast || 0;
  }
  B.prototype.get = function() {
    let g = {
      srcColorFactor: "SrcAlpha",
      dstColorFactor: "OneMinusSrcAlpha",
      func: "FuncAdd"
    };
    return this.additiveBlending && (g = {
      srcColorFactor: "SrcAlpha",
      dstColorFactor: "One",
      func: "FuncAdd"
    }), {
      blendCfg: g,
      opacity: this.opacity,
      color: {
        // Tonal corrections constants
        kGamma: this.kGamma,
        kSaturation: this.kSaturation,
        kBrightness: this.kBrightness,
        kContrast: this.kContrast,
        stretch: this.stretch,
        minCut: this.minCut,
        maxCut: this.maxCut,
        reversed: this.reversed,
        cmapName: this.colormap
      }
    };
  }, B.prototype.setBrightness = function(g) {
    g = +g || 0, this.kBrightness = Math.max(-1, Math.min(g, 1));
  }, B.prototype.getBrightness = function() {
    return this.kBrightness;
  }, B.prototype.setContrast = function(g) {
    g = +g || 0, this.kContrast = Math.max(-1, Math.min(g, 1));
  }, B.prototype.getContrast = function() {
    return this.kContrast;
  }, B.prototype.setSaturation = function(g) {
    g = +g || 0, this.kSaturation = Math.max(-1, Math.min(g, 1));
  }, B.prototype.getSaturation = function() {
    return this.kSaturation;
  }, B.prototype.setGamma = function(g) {
    g = +g, this.kGamma = Math.max(0.1, Math.min(g, 10));
  }, B.prototype.getGamma = function() {
    return this.kGamma;
  }, B.prototype.setOpacity = function(g) {
    g = +g, this.opacity = Math.max(0, Math.min(g, 1));
  }, B.prototype.setAlpha = B.prototype.setOpacity, B.prototype.getOpacity = function() {
    return this.opacity;
  }, B.prototype.getAlpha = B.prototype.getOpacity, B.prototype.setBlendingConfig = function(g = !1) {
    this.additiveBlending = g;
  }, B.prototype.getBlendingConfig = function() {
    return this.additiveBlending;
  };
  var A = function(g) {
    return g = g.toLowerCase(), B.COLORMAPS.includes(g) || (console.warn("The colormap '" + g + "' is not supported. You should use one of the following: " + B.COLORMAPS + `
'grayscale' has been chosen by default.`), g = "grayscale"), g;
  };
  return B.prototype.setColormap = function(g = "native", C = {}) {
    let I = A(g), E = C && C.stretch || this.stretch || "linear";
    E = E.toLowerCase();
    let o = !1;
    C && C.reversed === !0 && (o = !0), this.colormap = I, this.stretch = E, this.reversed = o;
  }, B.prototype.getColormap = function() {
    return this.colormap;
  }, B.prototype.setCuts = function(g, C) {
    this.minCut = g, this.maxCut = C;
  }, B.COLORMAPS = [], B;
}(), yQ = function() {
  let B = function(A, g) {
    this.id = "footprint-" + Z.uuidv4(), this.source = g, this.shapes = A, this.isShowing = !0;
  };
  return B.prototype.setCatalog = function(A) {
    this.source && this.source.setCatalog(A);
  }, B.prototype.show = function() {
    this.isShowing || (this.isShowing = !0, this.shapes.forEach((A) => A.show()));
  }, B.prototype.hide = function() {
    this.isShowing && (this.isShowing = !1, this.shapes.forEach((A) => A.hide()));
  }, B.prototype.select = function() {
    this.shapes.forEach((A) => A.select());
  }, B.prototype.deselect = function() {
    this.shapes.forEach((A) => A.deselect());
  }, B.prototype.setLineWidth = function(A) {
    this.shapes.forEach((g) => g.setLineWidth());
  }, B.prototype.setColor = function(A) {
    this.shapes.forEach((g) => g.setColor(A));
  }, B.prototype.setSelectionColor = function(A) {
    this.shapes.forEach((g) => g.setSelectionColor(A));
  }, B.prototype.isFootprint = function() {
    return !0;
  }, B.prototype.draw = function(A, g, C) {
    this.shapes.forEach((I) => I.draw(A, g, C));
  }, B.prototype.actionClicked = function() {
    this.source && this.source.actionClicked(), this.shapes.forEach((A) => A.select());
  }, B.prototype.actionOtherObjectClicked = function() {
    this.source && this.source.actionOtherObjectClicked(), this.shapes.forEach((A) => A.deselect());
  }, B.prototype.isInStroke = function(A, g, C, I) {
    return this.shapes.some((E) => E.isInStroke(A, g, C, I));
  }, B.prototype.getCatalog = function() {
    return this.source && this.source.catalog;
  }, B.prototype.intersectsBBox = function(A, g, C, I, E) {
    if (this.source) {
      let t = this.source;
      if (!t.isShowing)
        return !1;
      let s = null;
      if (t.x && t.y)
        s = {
          x: t.x,
          y: t.y
        };
      else {
        var o = jA.radecToViewXy(t.ra, t.dec, E);
        if (!o)
          return !1;
        s = {
          x: o[0],
          y: o[1]
        };
      }
      if (s.x >= A && s.x <= A + C && s.y >= g && s.y <= g + I)
        return !0;
    }
    return !1;
  }, B;
}(), eI = function() {
  function B(I, E, o, t, s) {
    this.aladin = I, this.options = I.options, this.aladinDiv = this.aladin.aladinDiv, this.popup = new dC(this.aladinDiv, this), this.createCanvases(), this.loadingState = !1;
    let M = this;
    try {
      const W = new NE(VA.wasmLibs.core, this.aladinDiv.id);
      this.aladin.wasm = W.webclient, this.wasm = this.aladin.wasm, CA.AL_USE_WASM.listenedBy(document.body, function(T) {
        let v = T.detail.callback;
        v(M.wasm);
      }), OI.COLORMAPS = this.wasm.getAvailableColormapList();
    } catch (W) {
      console.error(W), alert("Problem initializing Aladin Lite. Please contact the support by contacting Matthieu Baumann (baumannmatthieu0@gmail.com) or Thomas Boch (thomas.boch@astro.unistra.fr). You can also open an issue on the Aladin Lite github repository here: https://github.com/cds-astro/aladin-lite. Message error:" + W);
    }
    this.aladinDiv.ondrop = (W) => {
      Z.getDroppedFilesHandler(W).forEach((v) => {
        const _ = URL.createObjectURL(v);
        try {
          const a = M.aladin.createImageFITS(
            _,
            v.name,
            void 0,
            (iA, IA, z, oA) => {
              I.gotoRaDec(iA, IA), I.setFoV(z * 1.1);
            },
            void 0
          );
          M.setOverlayImageLayer(a, v.name);
        } catch (a) {
          let iA = zA.MOCFromURL(_);
          throw M.aladin.addMOC(iA), console.error("Only valid fits files supported (i.e. containig a WCS)", a), a;
        }
      });
    }, this.aladinDiv.ondragover = Z.dragOverHandler, this.location = E, this.fovDiv = o, this.mustClearCatalog = !0, this.mode = B.PAN, this.minFOV = this.maxFOV = null, this.healpixGrid = new HC(), this.then = Date.now();
    var N, k;
    N = k = 0, this.projection = uA.SIN, this.zoomFactor = this.wasm.getClipZoomFactor(), this.viewCenter = { lon: N, lat: k }, t ? this.cooFrame = t : this.cooFrame = RA.GAL, this.changeFrame(this.cooFrame);
    const U = 5e5, Y = 40;
    this.fovLimit = void 0;
    let K = s || 180;
    this.pinchZoomParameters = {
      isPinching: !1,
      // true if a pinch zoom is ongoing
      initialFov: void 0,
      initialDistance: void 0,
      initialAccDelta: Math.pow(U / K, 1 / Y)
    }, this.setZoom(K), this.imageLayers = /* @__PURE__ */ new Map(), this.overlayLayers = [], this.catalogs = [];
    var l = document.createElement("canvas");
    l.width = l.height = 24;
    var d = l.getContext("2d");
    d.lineWidth = 6, d.beginPath(), d.strokeStyle = "#eee", d.arc(12, 12, 8, 0, 2 * Math.PI, !0), d.stroke(), d.lineWidth = 3, d.beginPath(), d.strokeStyle = "#c38", d.arc(12, 12, 8, 0, 2 * Math.PI, !0), d.stroke(), this.catalogForPopup = zA.catalog({ shape: l, sourceSize: 24 }), this.catalogForPopup.hide(), this.catalogForPopup.setView(this), this.overlayForPopup = zA.graphicOverlay({ color: "#ee2345", lineWidth: 3 }), this.overlayForPopup.hide(), this.overlayForPopup.setView(this), this.overlays = [], this.mocs = [], this.allOverlayLayers = [], this.empty = !0, this.promises = [], this.firstHiPS = !0, this.curNorder = 1, this.realNorder = 1, this.imageLayersBeingQueried = /* @__PURE__ */ new Map(), this.dragging = !1, this.dragx = null, this.dragy = null, this.rightclickx = null, this.rightclicky = null, this.selectedLayer = "base", this.needRedraw = !0, this.fingersRotationParameters = {
      initialViewAngleFromCenter: void 0,
      initialFingerAngle: void 0,
      rotationInitiated: !1
    }, this.fadingLatestUpdate = null, this.dateRequestRedraw = null, C(this), this.resizeTimer = null;
    let q = new ResizeObserver(() => {
      M.fixLayoutDimensions(), M.requestRedraw();
    });
    this.throttledPositionChanged = Z.throttle(
      () => {
        var W = this.aladin.callbacksByEventName.positionChanged;
        if (typeof W == "function") {
          var T = this.aladin.pix2world(this.width / 2, this.height / 2);
          T !== void 0 && W({
            ra: T[0],
            dec: T[1],
            dragging: !0
          });
        }
      },
      B.CALLBACKS_THROTTLE_TIME_MS
    ), this.throttledZoomChanged = Z.throttle(
      () => {
        const W = this.fov;
        if (W !== this.oldFov) {
          const T = this.aladin.callbacksByEventName.zoomChanged;
          typeof T == "function" && T(W), this.oldFov = W;
        }
      },
      B.CALLBACKS_THROTTLE_TIME_MS
    ), q.observe(this.aladinDiv), M.fixLayoutDimensions(), M.requestRedraw();
  }
  B.PAN = 0, B.SELECT = 1, B.TOOL_SIMBAD_POINTER = 2, B.DRAW_SOURCES_WHILE_DRAGGING = !0, B.DRAW_MOCS_WHILE_DRAGGING = !0, B.CALLBACKS_THROTTLE_TIME_MS = 100, B.prototype.createCanvases = function() {
    var I = O(this.aladinDiv);
    I.find(".aladin-imageCanvas").remove(), I.find(".aladin-gridCanvas").remove(), I.find(".aladin-catalogCanvas").remove(), this.imageCanvas = O("<canvas class='aladin-imageCanvas'></canvas>").appendTo(this.aladinDiv)[0], this.gridCanvas = O("<canvas class='aladin-gridCanvas'></canvas>").appendTo(this.aladinDiv)[0], this.catalogCanvas = O("<canvas class='aladin-catalogCanvas'></canvas>").appendTo(this.aladinDiv)[0];
  }, B.prototype.fixLayoutDimensions = function() {
    Z.cssScale = void 0;
    var I = O(this.aladinDiv).width(), E = O(this.aladinDiv).height();
    this.width = Math.max(I, 1), this.height = Math.max(E, 1), this.cx = this.width / 2, this.cy = this.height / 2, this.largestDim = Math.max(this.width, this.height), this.smallestDim = Math.min(this.width, this.height), this.ratio = this.largestDim / this.smallestDim, this.mouseMoveIncrement = 160 / this.largestDim, this.imageCtx = this.imageCanvas.getContext("webgl2"), this.wasm.resize(this.width, this.height), this.catalogCtx = this.catalogCanvas.getContext("2d"), this.catalogCtx.canvas.width = this.width, this.catalogCtx.canvas.height = this.height, this.gridCtx = this.gridCanvas.getContext("2d"), this.gridCtx.canvas.width = this.width, this.gridCtx.canvas.height = this.height, A(this.imageCtx, this.aladin.options.pixelateCanvas), this.logoDiv || (this.logoDiv = O(this.aladinDiv).find(".aladin-logo")[0]), this.width > 800 ? (O(this.logoDiv).removeClass("aladin-logo-small"), O(this.logoDiv).addClass("aladin-logo-large"), O(this.logoDiv).css("width", "90px")) : (O(this.logoDiv).addClass("aladin-logo-small"), O(this.logoDiv).removeClass("aladin-logo-large"), O(this.logoDiv).css("width", "32px")), this.computeNorder(), this.redraw();
  };
  var A = function(I, E) {
    var o = !E;
    I.imageSmoothingEnabled = o, I.webkitImageSmoothingEnabled = o, I.mozImageSmoothingEnabled = o, I.msImageSmoothingEnabled = o, I.oImageSmoothingEnabled = o;
  };
  B.prototype.setMode = function(I) {
    this.mode = I, this.mode == B.SELECT ? this.setCursor("crosshair") : this.mode == B.TOOL_SIMBAD_POINTER ? (this.popup.hide(), this.catalogCanvas.style.cursor = "", O(this.catalogCanvas).addClass("aladin-sp-cursor")) : this.setCursor("default");
  }, B.prototype.setCursor = function(I) {
    this.catalogCanvas.style.cursor != I && this.mode != B.TOOL_SIMBAD_POINTER && (this.catalogCanvas.style.cursor = I);
  }, B.prototype.getCanvasDataURL = async function(I, E, o) {
    return function(s) {
      return new Promise((M, N) => {
        var k = new Image();
        k.src = s, k.onload = () => M(k), k.onerror = () => N(new Error("could not load image"));
      });
    }("data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAI4AAABTCAMAAAB+g8/LAAACx1BMVEVMaXGBYYSWk5i7ur1fGUW0Fzbi4OP////Qz9K2s7f////qyseffX7TxczMytBXU1ndrahOWXi0o7RaH0v///+1GjfYkY29srb///+1GTe0Fzajn6RgFkFdHkni3+GLV3PU0dXMubr6+vpmIktUJVKiGDqGcX7p5ujLwMJgFkFgFkFNOWnp1tZaHUi0FzaEZohkX2VVKVXUwcvy8vI4U4tQMWBXIk+NGT9ZIEx+Wn5vF0EUYqF3c3lgFkL5+PkUYqH///////9lFkG0FzYUYqFeNF/BwMs2WpP6+vrBv8JSJ1TNy85TJlO0FzaJhYsUYqF5GEEUYqF2Zo60FzazFza0FzYUYqGWdIsrWpWTGj6jGDp3Kk58Y4S0FzZgFkFXIU2OiY+vmqVhGENlGEJqQ2z///9SKFJTJlP///9pF0GOjpd0Ol6rFzi9sbm0Fza0FzYUYqGXmLp3TXJmHkhLSXy/jJBVK1ivrLDu7e7l5OYLCw6AYYRpFkGCIUYVYqGAZoqJfofez9hZPGtcW4phFkIUYqGVbG1BToTFw8ZqZGr4+PmIGkAWYqD6+vpaHUoUYqGEZoh5ZH2ceYAbGyCmFzmgGjsUYqGAYIOuiJJ3SW1PZJlNM0OliJ+MQF5uF0Gcmp8kXZpSKFWEZojDwcXq1tQzVY9pN2CyFzbZlZFHbKOZgpWjnaRlMlsUYqGHGD9FRElaHUiZfpfW1dddW2HMtsJ3k8NTJlPDT1WlMElcGkY6UYjMa2tDSH3IpKOEZoiFTWqni54DAwQsLDGsqa3Pu8cUFBnEtr8gHyU4Nz3cwsMKDA/GV1tGRUtCKjDczM7NfXzMvcza1Nv///+9PUmhfZRxY2y2KT/15eLo4ud5fKXCXmTnu7ekZ3pgFkFTJlOEZoiUGT5aHkp8GEBzF0G0FzadGDtKQnNeJ1JqJk5fGEReGkaDGT8UYqGlSw8iAAAAwXRSTlMA87vu8R/SwN6iQP7+/vf9/J75s4DT/v0gokr33vzj++7+9/Hz8/3u1tFw9P4f5nP9cvl0/vb+/vL79HH9++WPMFA7s1r++vRhscXEiWT9PvLQ+Ffzih/9/vb+9z3Enn7N/cWI/RDWPND+9/38gTx6uPj5/fn+/efauu7k8fnl0+ro/f33wvj7meDU2PeaZquWH9jJ1O0QrPfC0vXo+uHj+J7ETZvkpfzI+6e44qCorUr22cpX3xDd9VdUvtb6V9z+sGF5dwAACP1JREFUeF7s011r01AcBvATON8gFCgkV+2AFrKSm5MGCEKlDIqCgEgpXYUaOkanQLrtpupgCxTY8A3EDWToYBNlgFeiIOIX+f/z0pe96IcwSZtRxY0ByXaT3204nIfnPCHXLJFIJBKJgoe8LLyp/+fbPXJ16mvW3k7XsjiOs3xGd+1FoVAn12Hh1g7HqcYqMsdxGAZ0K8B15avOUkGPQymFvm0Plb6InrKOuqEbqoHVd1vPSfxk+fvT/VZRpBQ0aoLPtRW7VptRKD0VGTKcmNva/0biJPmVjDZUtXN8egKBXIM3IeC64NEohHlGvV6WxOcTj4hHhmq015dHyASh0ciXSKjUhAka5in21AMSi0ev3v7UEfEEjM5Rtbd+mPssSeQfz8JEIgZoR7VIHB6ubFvj4WqQ4xvnTqIkgE+j6KPQiSHOe54vlx0Krj38BYJ08bp27UUAcZyHQibiOJIsV9DXV4a1mrKYk8jFSndn+qCJwXuJZmYt2mKy6HvyemlJ8Zd7iSO3Bx8ANKCITDONQpTVtNCzam2vfHVBOK+OvLek/FRpmy4ABWBIob0X5TsF1Th6FY/NHC9NN5BOzadvzg5m06ldmGiSiQYAOCYwBpmNHyQaX+QW+ljbPDjkH5CJheCnnx+MDZU7j+FMcyqOSDU0Ye5jNL1UshhwaNvwo4SK4mYqNQjZGvzl/lkck1GKsPz7xiUu+0Nq2b+2VYVx/NDZJTYmnV2TpuvMsiJNhbSUZmMwSpssENJl7XSmrrDNpkpn3dqO4eraoqXFMmddBWcVncImDpgOMKiiImJu3t+Wl9a54UiccOxA8keY+5xzc25ugiTx+9s5fHL55D7nPM9dk5FY6NpO1wVgJ8g0pVIpv793mWLP31JEeiMKiCa5yeu8CRIeP8STySzLIMv5VSrl+e1YLne0Ap3BMMcnNE/XdV5Ybyer+lcOZyGeIsyKn+AxSDR8qcVwq9X6Lj+sDuwlm8FMJsiJ4o2fSX9fyeeXuY2D6MrpvDz1KEtylmIG/uh2Y6ZDlOomGxBaxx86CzovybniRG12VEEMUaCXLGV03svSPPaMXsBG8jKCDssHc3aE1BgLOj9OCzoshoYKdExxYL3zpTpuODZbo6+f7hKw0A5e5sBDqQ63MGcfwkxnHZXqeL+pQEd7kbpLdY5kwebt0f1HeGwbwYy8zsGMC7Ain9UfmE5va32pDqfXVuCjCwB73Vys0wUy+0f3fV6EeWLqkRn0U13QR9MTEOql4HXI5nZE304Ilo2E6KmkWnYCh9eKdMhI2LpxwU2xaYp10lZsdWKsbj138klVD/X55Q+Mnc/mOyC0bKLjvf3c4sBJB7mX8ekKdCb0rFpMh7ThrcPCNJhRK9kVrG/txkKGkMvHQe48wOpdu1dop6Q6j6N8Glxs8R9pgNAyXDSLdIJZyE4B+zkWS4QE7Fw33oyRYKxGyEWLYVTXmz/5jn+kGY0FRQYT8kp0tJPNfDb6AI6bpDrURtt/U6PRzArYTX5IaXZo+NzDGI+g99NE5/ivu5ebIbKxv1rEBhXpmL6F0yYn1YrqpDpjFHsHsCaKJUR9JwI66Dp5cY2fHaL3SZ75p3qd1QV4yLSDlkEr0mE2XcYQYF9RbHyzSMeaR66SpnS6GcmFrvzIVq2OthMgn9YyTP6cSawj2LhPJGCnrYAlxTrOeoROXSKH52umc2FfVTqsCFE9QgagAw6RztNuavNG8i7s5DE9wSIiHesuNNONP/ZKdFS5RXm1Oqtwo8KDhbGun0DIRXUKNlNGKab8HXRo8x5xYkyP8m1LQWcAVauj1QEz/AVC5jOkDHbk7mAzi9hsklr1ibAk04GBOksb4by2y8bRn1elw2rFqWACwLwOda6/WqTjXpnCyR6GGQAL7FWfuspuFk7aomRK9L+40lKzzhwUIQBNfzAOvOpgRqxzaOVvjCMi7HJc6N91gs7DE+M+OrWW9mSequ3tsFo19svymWwjFdlT0OF3dRGFIpkog1kEnZag0hfmSO4YX9u6UrOOqYcrSWic6LB4H5TDHENwdooSMB6/AfepNh2olTTpEh1jOUyJS3QCCU/uygCqUQfmeGmGz0p0wvfLYjGpTih9/ti1F1CtOvCVU5qwR/KZd7etLDbbIcHaz+euIVS7jiPAlYsKziiLr688tsSwhU877tu+XDyK/ofOxIZMHH3KD4m0D6q2QVpINu4p8lHyiQCRUCh6lYb2tUkZRJdI+5v+fCs38BGCyGgQaofHqC7DtrD4tx07aGkbDAM4/hTmB5gFhqAILAFs0SHYpqaMwkwRhtBWtmp0FobFURqw1uJlaQdO6SVMB0zZmNCeelLmbd1p32CXIjj2BNNkZUnyIZa0tKlujAFtveR3ed/b++fhvbwv/JcvDVFDmaSQg7YzSrkhile6MjW3OQQt4Ekkxp/PhsPJmRgDvZQp3mdlXVE4Bdo8tP36pqI0z/MP8d1T6FIdVWeXxEDW9TICPRUXfFwFzRzliZ0T/UnV63XqyhqL5Y77EXR58D5dW/KryUXXIfTY6TzBss2cNTsHdVlOIVIcRSPi3vq1lmNXdrx2guF548NbgJ4PR02lsG7mjEDHKCJP0/wen5hITEK3Y5crvY1oxRRC0HMHMyparudA1T0x0SmxTbqzaTTtzhvCaRx6blLwYTtnCv5paHPkbNSKGcuVDCF4BH1QXg50cuzx/GlzZO3iG5nO1jBcNIxCEPpjoyFhE0WSCgd/88IzZ/26kT++tq6MEItAv2yI2u4YoqZpiKR+8x+9ulB+TIiSTHKsjL+aVybGHEH/lEXMhRElUULUFZ1f94DlzfT0gntjJ5kVTX5JRZ0lKyclI8NAX00TGiKqhN9cUmSF06Mpmq7L2wHRxq5UFOXzyetMKA79RgQQ0TycCEgqpnRdJ/NsXkaU8kvnH4fvnSe9Oe9qfnXZ2I/DAHwq5cY0QrT4Ec0d4feLor5y8X14a+vycnExFotlQgwMSkQo+cRWD2EuLTve3LIh7L86fAaDFr/rbRgzXsuOz+fzFnNFo3AQZODWMJmCYdsPReDWMXEm2NTd4nA4HA6H4zc5mbo+QO8AVQAAAABJRU5ErkJggg==").then((s) => {
      I = I || "image/png";
      const M = this.wasm.canvas();
      var N = document.createElement("canvas");
      let k = window.devicePixelRatio;
      N.width = E || this.width * k, N.height = o || this.height * k;
      var U = N.getContext("2d");
      U.drawImage(M, 0, 0, N.width, N.height), U.drawImage(this.catalogCanvas, 0, 0, N.width, N.height);
      const Y = N.width - s.width, K = N.height - s.height;
      return U.drawImage(s, Y, K), N.toDataURL(I);
    });
  }, B.prototype.setActiveHiPSLayer = function(I) {
    if (!this.imageLayers.has(I))
      throw I + " does not exists. So cannot be selected";
    this.selectedLayer = I;
  };
  var g = function(I) {
    var E = !1;
    "ontouchstart" in window && (E = !0);
    let o = function(K) {
      const l = Z.relMouseCoords(I.imageCanvas, K);
      I.deselectObjects();
      try {
        const q = I.wasm.screenToWorld(l.x, l.y);
        var d = I.wasm.viewToICRSCooSys(q[0], q[1]);
        I.pointTo(d[0], d[1], { forceAnimation: !0 });
      } catch {
        return;
      }
    };
    E || O(I.catalogCanvas).dblclick(o), O(I.catalogCanvas).bind("contextmenu", function(K) {
      K.preventDefault();
    }, !1);
    let t = null, s = null;
    O(I.catalogCanvas).bind("mousedown touchstart", function(K) {
      K.preventDefault(), K.stopPropagation();
      const l = Z.relMouseCoords(I.imageCanvas, K);
      if (K.which === 3 || K.button === 2) {
        I.rightClick = !0, I.rightClickTimeStart = Date.now(), I.rightclickx = l.x, I.rightclicky = l.y;
        const q = I.imageLayers.get(I.selectedLayer);
        q && (q.imgFormat === "fits" ? (t = q.properties.minCutout || q.getColorCfg().minCut || 0, s = q.properties.maxCutout || q.getColorCfg().maxCut || 1) : (t = q.getColorCfg().minCut || 0, s = q.getColorCfg().maxCut || 1));
        return;
      }
      if (K.type === "touchstart" && K.originalEvent && K.originalEvent.targetTouches && K.originalEvent.targetTouches.length == 2) {
        I.dragging = !1, I.pinchZoomParameters.isPinching = !0;
        var d = I.wasm.getFieldOfView();
        I.pinchZoomParameters.initialFov = d, I.pinchZoomParameters.initialDistance = Math.sqrt(Math.pow(K.originalEvent.targetTouches[0].clientX - K.originalEvent.targetTouches[1].clientX, 2) + Math.pow(K.originalEvent.targetTouches[0].clientY - K.originalEvent.targetTouches[1].clientY, 2)), I.fingersRotationParameters.initialViewAngleFromCenter = I.wasm.getRotationAroundCenter(), I.fingersRotationParameters.initialFingerAngle = Math.atan2(K.originalEvent.targetTouches[1].clientY - K.originalEvent.targetTouches[0].clientY, K.originalEvent.targetTouches[1].clientX - K.originalEvent.targetTouches[0].clientX) * 180 / Math.PI;
        return;
      }
      return I.dragx = l.x, I.dragy = l.y, I.dragging = !0, I.mode == B.PAN ? I.setCursor("move") : I.mode == B.SELECT && (I.selectStartCoo = { x: I.dragx, y: I.dragy }), I.wasm.pressLeftMouseButton(I.dragx, I.dragy), !0;
    }), O(I.catalogCanvas).bind("mouseup", function(K) {
      if (I.rightClick) {
        Date.now() - I.rightClickTimeStart < 300 && I.aladin.contextMenu && I.aladin.contextMenu._showMenu(K), I.rightClick = !1, I.rightclickx = null, I.rightclicky = null, I.rightClickTimeStart = void 0;
        return;
      }
    }), O(I.catalogCanvas).bind("click mouseout touchend touchcancel", function(K) {
      if ((K.type === "touchend" || K.type === "touchcancel") && I.pinchZoomParameters.isPinching) {
        I.pinchZoomParameters.isPinching = !1, I.pinchZoomParameters.initialFov = I.pinchZoomParameters.initialDistance = void 0;
        return;
      }
      if ((K.type === "touchend" || K.type === "touchcancel") && I.fingersRotationParameters.rotationInitiated) {
        I.fingersRotationParameters.initialViewAngleFromCenter = void 0, I.fingersRotationParameters.initialFingerAngle = void 0, I.fingersRotationParameters.rotationInitiated = !1;
        return;
      }
      var l = I.realDragging === !0, d = I.mode === B.SELECT && I.dragging;
      if (I.dragging && (I.setCursor("default"), I.dragging = !1, l && (I.realDragging = !1)), d) {
        I.deselectObjects();
        const IA = I.getObjectsInBBox(
          I.selectStartCoo.x,
          I.selectStartCoo.y,
          I.dragx - I.selectStartCoo.x,
          I.dragy - I.selectStartCoo.y
        );
        if (IA.forEach((z) => {
          z.forEach((oA) => oA.select());
        }), IA.length > 0) {
          let z = {}, oA = IA.map((hA) => {
            let yA = hA[0].getCatalog(), rA = hA.map((nA) => nA instanceof yQ ? nA.source : nA), gA = {
              name: yA.name,
              color: yA.color,
              rows: rA,
              fields: yA.fields,
              fieldsClickedActions: yA.fieldsClickedActions
            };
            return yA.isObsCore && yA.isObsCore() && (z.save = !0), gA;
          });
          I.aladin.measurementTable.showMeasurement(oA, z);
        }
        I.selectedObjects = IA, I.aladin.fire(
          "selectend",
          IA
        ), I.requestRedraw();
        return;
      }
      I.mustClearCatalog = !0, I.dragx = I.dragy = null;
      const q = Z.relMouseCoords(I.imageCanvas, K);
      if ((K.type === "mouseout" || K.type === "touchend" || K.type === "touchcancel") && (I.updateLocation(q.x, q.y, !0), K.type === "mouseout")) {
        I.mode === B.TOOL_SIMBAD_POINTER && I.setMode(B.PAN);
        return;
      }
      if (I.mode == B.TOOL_SIMBAD_POINTER) {
        GQ(I, K);
        return;
      }
      var W = I.closestObjects(q.x, q.y, 5);
      if (!l && W) {
        I.deselectObjects();
        var T = W[0];
        T.marker ? (I.popup.setTitle(T.popupTitle), I.popup.setText(T.popupDesc), I.popup.setSource(T), I.popup.show()) : I.lastClickedObject && I.lastClickedObject.actionOtherObjectClicked && I.lastClickedObject.actionOtherObjectClicked(), T.actionClicked && T.actionClicked();
        var v = I.aladin.callbacksByEventName.objectClicked;
        if (typeof v == "function" && v(T), T.isFootprint()) {
          var _ = I.aladin.callbacksByEventName.footprintClicked;
          typeof _ == "function" && T != I.lastClickedObject && _(T);
        }
        I.lastClickedObject = T;
      } else if (!l && (I.deselectObjects(), I.lastClickedObject)) {
        I.aladin.measurementTable.hide(), I.popup.hide(), I.lastClickedObject instanceof MQ || I.lastClickedObject instanceof hQ || I.lastClickedObject instanceof eB ? I.lastClickedObject.deselect() : I.lastClickedObject.actionOtherObjectClicked();
        var v = I.aladin.callbacksByEventName.objectClicked;
        typeof v == "function" && v(null), I.lastClickedObject = null;
      }
      var a = I.aladin.callbacksByEventName.click;
      if (typeof a == "function") {
        var iA = I.aladin.pix2world(q.x, q.y);
        iA !== void 0 && a({ ra: iA[0], dec: iA[1], x: q.x, y: q.y, isDragging: l });
      }
      I.refreshProgressiveCats(), I.wasm.releaseLeftButtonMouse(q.x, q.y);
    });
    var M, N = null;
    O(I.catalogCanvas).bind("mousemove touchmove", function(K) {
      K.preventDefault();
      const l = Z.relMouseCoords(I.imageCanvas, K);
      if (I.rightClick) {
        var d = I.aladin.callbacksByEventName.rightClickMove;
        if (typeof d == "function") {
          d(l.x, l.y);
          return;
        }
        if (I.selectedLayer) {
          let rA = I.imageLayers.get(I.selectedLayer);
          const gA = {
            x: I.catalogCanvas.clientWidth * 0.5,
            y: I.catalogCanvas.clientHeight * 0.5
          }, nA = (l.x - gA.x) / I.catalogCanvas.clientWidth, og = -(l.y - gA.y) / I.catalogCanvas.clientHeight, PA = (s - t) * nA, gg = PA + (1 - 2 * og) * t, XA = PA + (1 + 2 * og) * s;
          gg <= XA && rA.setCuts(gg, XA);
        }
        return;
      }
      if (K.type === "touchmove" && I.pinchZoomParameters.isPinching && K.originalEvent && K.originalEvent.touches && K.originalEvent.touches.length == 2) {
        var q = Math.atan2(K.originalEvent.targetTouches[1].clientY - K.originalEvent.targetTouches[0].clientY, K.originalEvent.targetTouches[1].clientX - K.originalEvent.targetTouches[0].clientX) * 180 / Math.PI, W = I.fingersRotationParameters.initialFingerAngle - q;
        if (!I.fingersRotationParameters.rotationInitiated && Math.abs(W) >= 7 && (I.fingersRotationParameters.rotationInitiated = !0, I.fingersRotationParameters.initialFingerAngle = q, W = 0), I.fingersRotationParameters.rotationInitiated) {
          let nA = I.fingersRotationParameters.initialViewAngleFromCenter;
          I.wasm.getLongitudeReversed() ? nA -= W : nA += W, I.wasm.setRotationAroundCenter(nA);
        }
        const rA = Math.sqrt(Math.pow(K.originalEvent.touches[0].clientX - K.originalEvent.touches[1].clientX, 2) + Math.pow(K.originalEvent.touches[0].clientY - K.originalEvent.touches[1].clientY, 2)), gA = Math.min(Math.max(I.pinchZoomParameters.initialFov * I.pinchZoomParameters.initialDistance / rA, 2777777e-11), I.fovLimit);
        I.setZoom(gA);
        return;
      }
      if (!I.dragging && !I.moving && I.updateObjectsLookup(), (!I.dragging || E) && I.updateLocation(l.x, l.y, !1), !I.dragging) {
        var T = I.aladin.callbacksByEventName.mouseMove;
        if (typeof T == "function") {
          var v = I.aladin.pix2world(l.x, l.y);
          v !== void 0 ? T({ ra: v[0], dec: v[1], x: l.x, y: l.y }) : N != null && T({ ra: null, dec: null, x: l.x, y: l.y }), N = v;
        }
        if (!I.dragging && !I.mode == B.SELECT) {
          var _ = I.closestObjects(l.x, l.y, 5);
          if (_) {
            let rA = _[0];
            var a = I.aladin.callbacksByEventName.objectHovered, iA = I.aladin.callbacksByEventName.footprintHovered;
            I.setCursor("pointer"), typeof a == "function" && rA != M && a(rA), rA.isFootprint() && typeof iA == "function" && rA != M && iA(rA), M = rA;
          } else {
            I.setCursor("default");
            var IA = I.aladin.callbacksByEventName.objectHoveredStop;
            M && (M.isFootprint() && I.requestRedraw(), typeof IA == "function" && IA(M)), M = null;
          }
        }
        if (K.type === "mousemove")
          return;
      }
      if (!I.dragging)
        return;
      var z, oA;
      if (z = { x: I.dragx, y: I.dragy }, oA = { x: l.x, y: l.y }, I.dragx = l.x, I.dragy = l.y, I.mode == B.SELECT) {
        I.requestRedraw();
        return;
      }
      I.realDragging = !0, I.wasm.moveMouse(z.x, z.y, oA.x, oA.y), I.wasm.goFromTo(z.x, z.y, oA.x, oA.y);
      const [hA, yA] = I.wasm.getCenter();
      I.viewCenter.lon = hA, I.viewCenter.lat = yA, I.throttledPositionChanged();
    }), O(I.aladinDiv).onselectstart = function() {
      return !1;
    };
    var k = 0, U, Y;
    O(I.catalogCanvas).on("wheel", function(K) {
      if (K.preventDefault(), K.stopPropagation(), I.rightClick)
        return;
      var l = K.deltaY;
      K.hasOwnProperty("originalEvent") && (l = -K.originalEvent.deltaY);
      var d = Y || typeof Y < "u";
      d || (k === 0 && (U = (/* @__PURE__ */ new Date()).getTime()), k++, (/* @__PURE__ */ new Date()).getTime() - U > 100 && (k > 10 ? Y = !0 : Y = !1, d = !0));
      const q = (T) => {
        l > 0 ? I.increaseZoom(T) : I.decreaseZoom(T);
      };
      if (d && (Y ? ((/* @__PURE__ */ new Date()).getTime(), q(2e-3), (/* @__PURE__ */ new Date()).getTime()) : q(7e-3)), !I.debounceProgCatOnZoom) {
        var W = I;
        I.debounceProgCatOnZoom = Z.debounce(function() {
          W.refreshProgressiveCats(), W.drawAllOverlays();
        }, 300);
      }
      return I.debounceProgCatOnZoom(), I.throttledZoomChanged(), !1;
    });
  }, C = function(I) {
    var E = new xC();
    E.domElement.style.top = "50px", O("#aladin-statsDiv").length > 0 && O("#aladin-statsDiv")[0].appendChild(E.domElement), I.stats = E, g(I), I.displayHpxGrid = !1, I.displayCatalog = !1, I.displayReticle = !0;
  };
  return B.prototype.updateLocation = function(I, E, o) {
    if (o)
      this.location.update(this.viewCenter.lon, this.viewCenter.lat, this.cooFrame, !0);
    else {
      let t = this.wasm.screenToWorld(I, E);
      t && (t[0] < 0 && (t = [t[0] + 360, t[1]]), this.location.update(t[0], t[1], this.cooFrame, !1));
    }
  }, B.prototype.requestRedrawAtDate = function(I) {
    this.dateRequestDraw = I;
  }, B.prototype.getViewParams = function() {
    var I = this.width > this.height ? this.fov / this.width : this.fov / this.height;
    return {
      fov: [this.width * I, this.height * I],
      width: this.width,
      height: this.height
    };
  }, B.FPS_INTERVAL = 1e3 / 140, B.prototype.redraw = function() {
    const I = Date.now(), E = I - this.then;
    try {
      this.moving = this.wasm.update(E);
    } catch (t) {
      console.warn(t);
    }
    (this.wasm.isRendering() || this.needRedraw) && this.drawAllOverlays(), this.needRedraw = !1, this.then = I, mC(this.redraw.bind(this));
  }, B.prototype.drawAllOverlays = function() {
    var I = this.catalogCtx, E = !1;
    if (this.mustClearCatalog && (I.clearRect(0, 0, this.width, this.height), E = !0, this.mustClearCatalog = !1), this.catalogs && this.catalogs.length > 0 && this.displayCatalog && (!this.dragging || B.DRAW_SOURCES_WHILE_DRAGGING)) {
      E || (I.clearRect(0, 0, this.width, this.height), E = !0);
      for (var o = 0; o < this.catalogs.length; o++) {
        var t = this.catalogs[o];
        t.draw(I, this.cooFrame, this.width, this.height, this.largestDim, this.zoomFactor);
      }
    }
    this.catalogForPopup.isShowing && this.catalogForPopup.sources.length > 0 && (E || (I.clearRect(0, 0, this.width, this.height), E = !0), this.catalogForPopup.draw(I, this.cooFrame, this.width, this.height, this.largestDim, this.zoomFactor), this.overlayForPopup.isShowing && this.overlayForPopup.draw(I, this.cooFrame, this.width, this.height, this.largestDim, this.zoomFactor));
    var s = this.catalogCtx;
    if (this.overlays && this.overlays.length > 0 && (!this.dragging || B.DRAW_SOURCES_WHILE_DRAGGING)) {
      E || (I.clearRect(0, 0, this.width, this.height), E = !0);
      for (var o = 0; o < this.overlays.length; o++)
        this.overlays[o].draw(s);
    }
    var M = I;
    if (this.displayHpxGrid) {
      E || (I.clearRect(0, 0, this.width, this.height), E = !0);
      var N = this.getVisibleCells(3), k = null;
      this.curNorder >= 3 && (this.curNorder == 3 ? k = N : k = this.getVisibleCells(this.curNorder)), k && this.curNorder > 3 ? this.healpixGrid.redraw(M, k, this.fov, this.curNorder) : this.healpixGrid.redraw(M, N, this.fov, 3);
    }
    var U = I;
    if (this.mode == B.SELECT) {
      if (this.dragging) {
        E || (U.clearRect(0, 0, this.width, this.height), E = !0), U.fillStyle = "rgba(100, 240, 110, 0.25)";
        var Y = this.dragx - this.selectStartCoo.x, K = this.dragy - this.selectStartCoo.y;
        U.fillRect(this.selectStartCoo.x, this.selectStartCoo.y, Y, K);
      }
    } else if (this.displayReticle) {
      if (E || (I.clearRect(0, 0, this.width, this.height), E = !0), !this.reticleCache) {
        var l = document.createElement("canvas"), d = this.options.reticleSize;
        l.width = d, l.height = d;
        var q = l.getContext("2d");
        q.lineWidth = 2, q.strokeStyle = this.options.reticleColor, q.beginPath(), q.moveTo(d / 2, d / 2 + (d / 2 - 1)), q.lineTo(d / 2, d / 2 + 2), q.moveTo(d / 2, d / 2 - (d / 2 - 1)), q.lineTo(d / 2, d / 2 - 2), q.moveTo(d / 2 + (d / 2 - 1), d / 2), q.lineTo(d / 2 + 2, d / 2), q.moveTo(d / 2 - (d / 2 - 1), d / 2), q.lineTo(d / 2 - 2, d / 2), q.stroke(), this.reticleCache = l;
      }
      U.drawImage(this.reticleCache, this.width / 2 - this.reticleCache.width / 2, this.height / 2 - this.reticleCache.height / 2);
    }
    if (this.projection == uA.SIN && this.fov >= 60 && this.aladin.options.showAllskyRing === !0) {
      E || (U.clearRect(0, 0, this.width, this.height), E = !0), U.strokeStyle = this.aladin.options.allskyRingColor;
      var W = this.aladin.options.allskyRingWidth;
      U.lineWidth = W, U.beginPath();
      const v = (this.cy - W / 2 + 1) / this.zoomFactor;
      U.arc(this.cx, this.cy, v, 0, 2 * Math.PI, !0), U.stroke();
    }
  }, B.prototype.refreshProgressiveCats = function() {
    if (this.catalogs)
      for (var I = 0; I < this.catalogs.length; I++)
        this.catalogs[I].type == "progressivecat" && this.catalogs[I].loadNeededTiles();
  }, B.prototype.getVisiblePixList = function(I) {
    var E = [];
    let o = this.wasm.screenToWorld(this.cx, this.cy);
    const [t, s] = this.wasm.viewToICRSCooSys(o[0], o[1]);
    var M = this.fov * 0.5 * this.ratio;
    return this.wasm.queryDisc(I, t, s, M).forEach((N) => E.push(Number(N))), E;
  }, B.prototype.deselectObjects = function() {
    this.selectedObjects && (this.selectedObjects.forEach((I) => {
      I.forEach((E) => E.deselect());
    }), this.aladin.measurementTable.hide(), this.selectedObjects = null);
  }, B.prototype.getVisibleCells = function(I) {
    return this.wasm.getVisibleCells(I);
  }, B.prototype.setZoom = function(I) {
    this.wasm.setFieldOfView(I), this.updateZoomState();
  }, B.prototype.increaseZoom = function(I) {
    let t = this.pinchZoomParameters.initialAccDelta + I, s = 5e5 / Math.pow(t, 40);
    s < 2777777e-11 && (s = 2777777e-11), this.pinchZoomParameters.initialAccDelta = t, this.setZoom(s);
  }, B.prototype.decreaseZoom = function(I) {
    let t = this.pinchZoomParameters.initialAccDelta - I;
    t <= 0 && (t = 1e-3);
    let s = 5e5 / Math.pow(t, 40);
    s >= this.fovLimit && (s = this.fovLimit), this.pinchZoomParameters.initialAccDelta = t, this.setZoom(s);
  }, B.prototype.setRotation = function(I) {
    this.wasm.setRotationAroundCenter(I);
  }, B.prototype.setGridConfig = function(I) {
    this.wasm.setGridConfig(I), I && (I.hasOwnProperty("enabled") && (I.enabled === !0 ? CA.COO_GRID_ENABLED.dispatchedTo(this.aladinDiv) : CA.COO_GRID_DISABLED.dispatchedTo(this.aladinDiv)), I.color && CA.COO_GRID_UPDATED.dispatchedTo(this.aladinDiv, { color: I.color, opacity: I.opacity })), this.requestRedraw();
  }, B.prototype.updateZoomState = function() {
    this.zoomFactor = this.wasm.getClipZoomFactor();
    let I = this.wasm.getFieldOfView();
    const E = 5e5, o = 40;
    if (this.pinchZoomParameters.initialAccDelta = Math.pow(E / I, 1 / o), this.fov = I, this.computeNorder(), isNaN(this.fov)) {
      this.fovDiv.html("FoV:");
      return;
    }
    var t;
    this.projection.fov <= I && (I = this.projection.fov), I > 1 ? t = Math.round(I * 100) / 100 + "°" : I * 60 > 1 ? t = Math.round(I * 60 * 100) / 100 + "'" : t = Math.round(I * 3600 * 100) / 100 + '"', this.fovDiv.html("FoV: " + t);
  }, B.prototype.computeNorder = function() {
    var I = this.wasm.getNOrder();
    this.realNorder = I, this.fov <= 50 && I <= 2 && (I = 3), this.curNorder = I;
  }, B.prototype.untaintCanvases = function() {
    this.createCanvases(), g(this), this.fixLayoutDimensions();
  }, B.prototype.setOverlayImageLayer = function(I, E = "overlay") {
    return this.imageLayersBeingQueried.set(E, I), this.addImageLayer(I, E), I;
  }, B.prototype.addLayer = function(I) {
    const E = I.layer;
    this.overlayLayers.findIndex((M) => M == E) == -1 && this.overlayLayers.push(E), this.options.log && cB.log("setImageLayer", I.url);
    const t = this.overlayLayers[this.overlayLayers.length - 1];
    this.selectedLayer = t;
    let s = this.imageLayers.get(E);
    s && (s.added = !1), this.imageLayers.set(E, I), CA.HIPS_LAYER_ADDED.dispatchedTo(this.aladinDiv, { layer: I });
  }, B.prototype.addImageLayer = function(I, E) {
    let o = this;
    const t = I.query;
    this.promises.push(t), Promise.allSettled(this.promises).then(() => t).then((s) => {
      const M = s.add(E);
      return o.loadingState = !0, CA.LOADING_STATE.dispatchedTo(this.aladinDiv, { loading: !0 }), M;
    }).then((s) => {
      this.empty = !1, s.children ? s.children.forEach((M) => {
        this.addLayer(M);
      }) : this.addLayer(s);
    }).catch((s) => {
      throw s;
    }).finally(() => {
      o.loadingState = !1, CA.LOADING_STATE.dispatchedTo(this.aladinDiv, { loading: !1 }), o.imageLayersBeingQueried.delete(E);
      let s = this.promises.findIndex((N) => N == t);
      if (this.promises.splice(s, 1), this.promises.length === 0)
        if (o.empty) {
          const N = Math.round(Math.random()), k = VA.DEFAULT_OPTIONS.surveyUrl[N];
          o.aladin.setBaseImageLayer(k);
        } else
          o.renameLayer(this.overlayLayers[0], "base");
    });
  }, B.prototype.renameLayer = function(I, E) {
    if (I === E)
      return;
    this.wasm.renameLayer(I, E);
    let o = this.imageLayers.get(I);
    o.layer = E;
    const t = this.overlayLayers.findIndex((s) => s == I);
    this.overlayLayers[t] = E, this.imageLayers.delete(I), this.imageLayers.set(E, o), this.selectedLayer === I && (this.selectedLayer = E), CA.HIPS_LAYER_RENAMED.dispatchedTo(this.aladinDiv, { layer: I, newLayer: E });
  }, B.prototype.swapLayers = function(I, E) {
    this.wasm.swapLayers(I, E);
    const o = this.overlayLayers.findIndex((M) => M == I), t = this.overlayLayers.findIndex((M) => M == E), s = this.overlayLayers[o];
    this.overlayLayers[o] = this.overlayLayers[t], this.overlayLayers[t] = s, CA.HIPS_LAYER_SWAP.dispatchedTo(this.aladinDiv, { firstLayer: I, secondLayer: E });
  }, B.prototype.removeImageLayer = function(I) {
    let E = this.imageLayers.get(I);
    if (E === void 0)
      return;
    E.added && this.wasm.removeLayer(I), E.added = !1;
    const o = this.overlayLayers.findIndex((s) => s == I);
    if (o == -1)
      return;
    if (this.imageLayers.delete(I), this.overlayLayers.splice(o, 1), this.overlayLayers.length === 0)
      this.empty = !0, this.selectedLayer = "base";
    else if (this.selectedLayer === I) {
      const s = this.overlayLayers[this.overlayLayers.length - 1];
      this.selectedLayer = s;
    }
    if (CA.HIPS_LAYER_REMOVED.dispatchedTo(this.aladinDiv, { layer: I }), this.promises.length === 0 && this.empty) {
      const s = Math.round(Math.random()), M = VA.DEFAULT_OPTIONS.surveyUrl[s];
      this.aladin.setBaseImageLayer(M);
    }
  }, B.prototype.setHiPSUrl = function(I, E) {
    try {
      this.wasm.setHiPSUrl(I, E);
    } catch (o) {
      console.error(o);
    }
  }, B.prototype.getImageLayer = function(I = "base") {
    let E = this.imageLayersBeingQueried.get(I), o = this.imageLayers.get(I);
    return E || o;
  }, B.prototype.requestRedraw = function() {
    this.needRedraw = !0;
  }, B.prototype.setProjection = function(I) {
    switch (this.fovLimit = 1e3, I) {
      case "TAN":
        this.projection = uA.TAN, this.fovLimit = 180;
        break;
      case "STG":
        this.projection = uA.STG;
        break;
      case "SIN":
        this.projection = uA.SIN;
        break;
      case "ZEA":
        this.projection = uA.ZEA;
        break;
      case "FEYE":
        this.projection = uA.FEYE;
        break;
      case "AIR":
        this.projection = uA.AIR;
        break;
      case "ARC":
        this.projection = uA.ARC;
        break;
      case "NCP":
        this.projection = uA.NCP;
        break;
      case "ARC":
        this.projection = uA.ARC;
        break;
      case "ZEA":
        this.projection = uA.ZEA;
        break;
      case "MOL":
        this.projection = uA.MOL;
        break;
      case "AIT":
        this.projection = uA.AIT;
        break;
      case "MER":
        this.projection = uA.MER, this.fovLimit = 360;
        break;
      case "CAR":
        this.projection = uA.CAR, this.fovLimit = 360;
        break;
      case "CEA":
        this.projection = uA.CEA, this.fovLimit = 360;
        break;
      case "CYP":
        this.projection = uA.CYP, this.fovLimit = 360;
        break;
      case "COD":
        this.projection = uA.COD;
        break;
      case "HPX":
        this.projection = uA.HPX, this.fovLimit = 360;
        break;
    }
    this.wasm.setProjection(I), this.updateZoomState(), this.requestRedraw();
  }, B.prototype.changeFrame = function(I) {
    this.cooFrame = I, this.cooFrame.system == RA.SYSTEMS.GAL ? this.wasm.setCooSystem(VA.wasmLibs.core.CooSystem.GAL) : this.cooFrame.system == RA.SYSTEMS.J2000 && this.wasm.setCooSystem(VA.wasmLibs.core.CooSystem.ICRS), this.cooFrame.label == "J2000d" ? this.setGridConfig({ fmt: "HMS" }) : this.setGridConfig({ fmt: "DMS" });
    let [E, o] = this.wasm.getCenter();
    this.viewCenter.lon = E, this.viewCenter.lat = o, this.location.update(this.viewCenter.lon, this.viewCenter.lat, this.cooFrame, !0), this.requestRedraw();
  }, B.prototype.showHealpixGrid = function(I) {
    this.displayHpxGrid = I, this.displayHpxGrid || (this.mustClearCatalog = !0), this.requestRedraw();
  }, B.prototype.showSurvey = function(I) {
    this.getImageLayer().setAlpha(I ? 1 : 0), this.requestRedraw();
  }, B.prototype.showCatalog = function(I) {
    this.displayCatalog = I, this.displayCatalog || (this.mustClearCatalog = !0), this.requestRedraw();
  }, B.prototype.showReticle = function(I) {
    this.displayReticle = I, this.displayReticle || (this.mustClearCatalog = !0), this.requestRedraw();
  }, B.prototype.pointTo = function(I, E, o) {
    if (I = parseFloat(I), E = parseFloat(E), !(isNaN(I) || isNaN(E))) {
      this.viewCenter.lon = I, this.viewCenter.lat = E, this.location.update(this.viewCenter.lon, this.viewCenter.lat, this.cooFrame, !0), this.wasm.setCenter(this.viewCenter.lon, this.viewCenter.lat), this.requestRedraw();
      var t = this;
      setTimeout(function() {
        t.refreshProgressiveCats();
      }, 1e3);
    }
  }, B.prototype.makeUniqLayerName = function(I) {
    if (!this.layerNameExists(I))
      return I;
    for (var E = 1; ; ++E) {
      var o = I + "_" + E;
      if (!this.layerNameExists(o))
        return o;
    }
  }, B.prototype.layerNameExists = function(I) {
    for (var E = this.allOverlayLayers, o = 0; o < E.length; o++)
      if (I == E[o].name)
        return !0;
    return !1;
  }, B.prototype.removeLayers = function() {
    this.catalogs = [], this.overlays = [], this.mocs = [], this.allOverlayLayers = [], this.requestRedraw();
  }, B.prototype.removeLayer = function(I) {
    let E = this.allOverlayLayers.indexOf(I);
    this.allOverlayLayers.splice(E, 1), I.type == "catalog" || I.type == "progressivecat" ? (E = this.catalogs.indexOf(I), this.catalogs.splice(E, 1)) : I.type == "moc" ? (E = this.mocs.indexOf(I), this.mocs.splice(E, 1)[0].delete()) : I.type == "overlay" && (E = this.overlays.indexOf(I), this.overlays.splice(E, 1)), CA.GRAPHIC_OVERLAY_LAYER_REMOVED.dispatchedTo(this.aladinDiv, { layer: I }), this.requestRedraw();
  }, B.prototype.addCatalog = function(I) {
    I.name = this.makeUniqLayerName(I.name), this.allOverlayLayers.push(I), this.catalogs.push(I), I.type == "catalog" ? I.setView(this) : I.type == "progressivecat" && I.init(this);
  }, B.prototype.addOverlay = function(I) {
    I.name = this.makeUniqLayerName(I.name), this.overlays.push(I), this.allOverlayLayers.push(I), I.setView(this);
  }, B.prototype.addMOC = function(I) {
    I.name = this.makeUniqLayerName(I.name), I.setView(this);
  }, B.prototype.getObjectsInBBox = function(I, E, o, t) {
    o < 0 && (I = I + o, o = -o), t < 0 && (E = E + t, t = -t);
    var s = [], M, N, k, U, Y, K = [];
    if (this.catalogs) {
      for (var l = 0; l < this.catalogs.length; l++)
        if (M = this.catalogs[l], !!M.isShowing) {
          N = M.getSources();
          for (var d = 0; d < N.length; d++)
            k = N[d], !(!k.isShowing || !k.x || !k.y) && k.x >= I && k.x <= I + o && k.y >= E && k.y <= E + t && K.push(k);
          if (U = M.getFootprints(), U)
            for (var d = 0; d < U.length; d++)
              Y = U[d], Y.intersectsBBox(I, E, o, t, this) && K.push(Y);
          K.length > 0 && s.push(K), K = [];
        }
    }
    return s;
  }, B.prototype.updateObjectsLookup = function() {
    this.objLookup = [];
    var I, E, o, t, s;
    if (this.catalogs) {
      for (var M = 0; M < this.catalogs.length; M++)
        if (I = this.catalogs[M], !!I.isShowing) {
          E = I.getSources();
          for (var N = 0; N < E.length; N++)
            o = E[N], !(!o.isShowing || !o.x || !o.y) && (t = Math.round(o.x), s = Math.round(o.y), typeof this.objLookup[t] > "u" && (this.objLookup[t] = []), typeof this.objLookup[t][s] > "u" && (this.objLookup[t][s] = []), this.objLookup[t][s].push(o));
        }
    }
  }, B.prototype.closestFootprints = function(I, E, o, t) {
    if (!I)
      return null;
    let s = null;
    return I.forEach((M) => {
      if (M.isShowing && M.isInStroke(E, this, o, t)) {
        s = M;
        return;
      }
    }), s;
  }, B.prototype.closestObjects = function(I, E, o) {
    var t, s = this.catalogCanvas, M = s.getContext("2d");
    let N = M.lineWidth;
    if (M.lineWidth = 6, this.overlays)
      for (var k = 0; k < this.overlays.length; k++) {
        t = this.overlays[k];
        let W = this.closestFootprints(t.overlayItems, M, I, E);
        if (W)
          return M.lineWidth = N, [W];
      }
    if (this.catalogs)
      for (var k = 0; k < this.catalogs.length; k++) {
        let T = this.catalogs[k], v = this.closestFootprints(T.footprints, M, I, E);
        if (v)
          return M.lineWidth = N, [v];
      }
    if (!this.objLookup)
      return M.lineWidth = N, null;
    M.lineWidth = N;
    for (var U, Y, K = 0; K <= o; K++) {
      U = Y = null;
      for (var l = -o; l <= o; l++)
        if (this.objLookup[I + l]) {
          for (var d = -o; d <= o; d++)
            if (this.objLookup[I + l][E + d]) {
              var q = l * l + d * d;
              (!U || q < Y) && (U = this.objLookup[I + l][E + d], Y = q);
            }
        }
      if (U)
        return U;
    }
    return null;
  }, B;
}(), GB = function() {
  let B = function(A, g, C, I) {
    this.ra = A, this.dec = g, this.data = C, this.catalog = null, this.marker = I && I.marker || !1, this.marker && (this.popupTitle = I && I.popupTitle ? I.popupTitle : "", this.popupDesc = I && I.popupDesc ? I.popupDesc : "", this.useMarkerDefaultIcon = I && I.useMarkerDefaultIcon !== void 0 ? I.useMarkerDefaultIcon : !0), this.isShowing = !0, this.isSelected = !1;
  };
  return B.prototype.setCatalog = function(A) {
    this.catalog = A;
  }, B.prototype.getCatalog = function() {
    return this.catalog;
  }, B.prototype.show = function() {
    this.isShowing || (this.isShowing = !0, this.catalog && this.catalog.reportChange());
  }, B.prototype.hide = function() {
    this.isShowing && (this.isShowing = !1, this.catalog && this.catalog.reportChange());
  }, B.prototype.select = function() {
    this.isSelected || (this.isSelected = !0, this.catalog && this.catalog.reportChange());
  }, B.prototype.deselect = function() {
    this.isSelected && (this.isSelected = !1, this.catalog && this.catalog.reportChange());
  }, B.prototype.actionClicked = function() {
    if (this.catalog && this.catalog.onClick) {
      var A = this.catalog.view;
      if (this.catalog.onClick == "showTable") {
        this.select();
        let I = {
          rows: [this],
          fields: this.catalog.fields,
          fieldsClickedActions: this.catalog.fieldsClickedActions,
          name: this.catalog.name,
          color: this.catalog.color
        }, E = {};
        E.save = !0, A.aladin.measurementTable.hide(), A.aladin.measurementTable.showMeasurement([I], E);
      } else if (this.catalog.onClick == "showPopup") {
        A.popup.setTitle("<br><br>");
        var g = '<div class="aladin-marker-measurement">';
        g += "<table>";
        for (var C in this.data)
          g += "<tr><td>" + C + "</td><td>" + this.data[C] + "</td></tr>";
        g += "</table>", g += "</div>", A.popup.setText(g), A.popup.setSource(this), A.popup.show();
      } else
        typeof this.catalog.onClick == "function" && (this.catalog.onClick(this), A.lastClickedObject = this);
    }
  }, B.prototype.isFootprint = function() {
    return !1;
  }, B.prototype.actionOtherObjectClicked = function() {
    this.catalog && this.catalog.onClick && this.deselect();
  }, B;
}(), _A = function() {
  let B = {};
  return B.curIdx = 0, B.colors = ["#ff0000", "#0000ff", "#99cc00", "#ffff00", "#000066", "#00ffff", "#9900cc", "#0099cc", "#cc9900", "#cc0099", "#00cc99", "#663333", "#ffcc9a", "#ff9acc", "#ccff33", "#660000", "#ffcc33", "#ff00ff", "#00ff00", "#ffffff"], B.standardizedColors = {}, B.getNextColor = function() {
    var A = B.colors[B.curIdx % B.colors.length];
    return B.curIdx++, A;
  }, B.getLabelColorForBackground = function(A) {
    var g = "#eee", C = "#111", I = A.match(/^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/);
    if (I == null)
      return C;
    var E = parseInt(I[1]), o = parseInt(I[2]), t = parseInt(I[3]), s = 1 - (0.299 * E + 0.587 * o + 0.114 * t) / 255;
    return s < 0.5 ? C : g;
  }, B.hexToRgb = function(A) {
    var g = /^#?([a-f\d])([a-f\d])([a-f\d])$/i;
    A = A.replace(g, function(I, E, o, t) {
      return E + E + o + o + t + t;
    });
    var C = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(A);
    return C ? {
      r: parseInt(C[1], 16),
      g: parseInt(C[2], 16),
      b: parseInt(C[3], 16)
    } : null;
  }, B.hexToRgba = function(A) {
    var g = /^#?([a-f\d])([a-f\d])([a-f\d])([a-f\d])$/i;
    A = A.replace(g, function(I, E, o, t, s) {
      return E + E + o + o + t + t + s + s;
    });
    var C = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(A);
    return C ? {
      r: parseInt(C[1], 16),
      g: parseInt(C[2], 16),
      b: parseInt(C[3], 16),
      a: parseInt(C[4], 16)
    } : null;
  }, B.rgbToHex = function(A, g, C) {
    return "#" + ((1 << 24) + (A << 16) + (g << 8) + C).toString(16).slice(1);
  }, B.standardizeColor = function(A) {
    if (A.match(/^#?[0-9A-Fa-f]{6}$/) || A.match(/^#?[0-9A-Fa-f]{3}$/))
      return A;
    if (!B.standardizedColors[A]) {
      var g = document.createElement("canvas").getContext("2d");
      g.fillStyle = A;
      const C = g.fillStyle;
      B.standardizedColors[A] = C;
    }
    return B.standardizedColors[A];
  }, B;
}(), kE = function() {
  function B() {
  }
  return B.handleActions = function(A, g) {
    hB.parse(
      A,
      (C, I) => {
        let E = [];
        I.forEach((t) => {
          let s = {};
          for (const [N, k] of Object.entries(C)) {
            var M = k.name;
            s[M] = t[k.idx];
          }
          E.push({ data: s });
        });
        let o = {
          name: "Datalink:" + A,
          color: "purple",
          rows: E,
          fields: C,
          fieldsClickedActions: {
            access_url: (t) => {
              let s = t.access_url, M = t.content_type, N = t.content_qualifier, k = () => {
                var U = (K, l, d, q) => {
                  g.gotoRaDec(K, l), g.setFoV(d);
                };
                let Y = g.createImageFITS(s, s, {}, U);
                g.setOverlayImageLayer(Y, Z.uuidv4());
              };
              switch (M) {
                case "application/hips":
                  let U = g.newImageSurvey(s);
                  g.setOverlayImageLayer(U, Z.uuidv4());
                  break;
                case "application/fits":
                  N === "cube" && console.warn("Cube not handled, only first slice downloaded"), k();
                  break;
                case "image/fits":
                  N === "cube" && console.warn("Cube not handled, only first slice downloaded"), k();
                  break;
                default:
                  Z.download(s);
                  break;
              }
            }
          }
        };
        g.measurementTable.showMeasurement([o], { save: !0 });
      }
    );
  }, B;
}(), rQ = function() {
  B.MANDATORY_FIELDS = {
    dataproduct_type: { name: "dataproduct_type", ucd: "meta.id", utype: "ObsDataset.dataProductType", units: null },
    calib_level: { name: "calib_level", ucd: "meta.code;obs.calib", utype: "ObsDataset.calibLevel", units: null },
    obs_collection: { name: "obs_collection", ucd: "meta.id", utype: "DataID.collection", units: null },
    obs_id: { name: "obs_id", ucd: "meta.id", utype: "DataID.observationID", units: null },
    obs_publisher_did: { name: "obs_publisher_did", ucd: "meta.ref.uri;meta.curation", utype: "Curation.publisherDID", units: null },
    access_url: { name: "access_url", ucd: "meta.ref.url", utype: "Access.reference", units: null },
    access_format: { name: "access_format", ucd: "meta.code.mime", utype: "Access.format", units: null },
    access_estsize: { name: "access_estsize", ucd: "phys.size;meta.file", utype: "Access.size", units: "kbyte" },
    target_name: { name: "target_name", ucd: "meta.id;src", utype: "Target.name", units: null },
    s_ra: { name: "s_ra", ucd: "pos.eq.ra", utype: "Char.SpatialAxis.Coverage.Location.Coord.Position2D.Value2.C1", units: "deg" },
    s_dec: { name: "s_dec", ucd: "pos.eq.dec", utype: "Char.SpatialAxis.Coverage.Location.Coord.Position2D.Value2.C2", units: "deg" },
    s_fov: { name: "s_fov", ucd: "phys.angSize;instr.fov", utype: "Char.SpatialAxis.Coverage.Bounds.Extent.diameter", units: "deg" },
    s_region: { name: "s_region", ucd: "pos.outline;obs.field", utype: "Char.SpatialAxis.Coverage.Support.Area", units: null },
    s_resolution: { name: "s_resolution", ucd: "pos.angResolution", utype: "Char.SpatialAxis.Resolution.Refval.value", units: "arcsec" },
    s_xel1: { name: "s_xel1", ucd: "meta.number", utype: "Char.SpatialAxis.numBins1", units: null },
    s_xel2: { name: "s_xel2", ucd: "meta.number", utype: "Char.SpatialAxis.numBins2", units: null },
    t_min: { name: "t_min", ucd: "time.start;obs.exposure", utype: "Char.TimeAxis.Coverage.Bounds.Limits.StartTime", units: "d" },
    t_max: { name: "t_max", ucd: "time.end;obs.exposure", utype: "Char.TimeAxis.Coverage.Bounds.Limits.StopTime", units: "d" },
    t_exptime: { name: "t_exptime", ucd: "time.duration;obs.exposure", utype: "Char.TimeAxis.Coverage.Support.Extent", units: "s" },
    t_resolution: { name: "t_resolution", ucd: "time.resolution", utype: "Char.TimeAxis.Resolution.Refval.value", units: "s" },
    t_xel: { name: "t_xel", ucd: "meta.number", utype: "Char.TimeAxis.numBins", units: null },
    em_min: { name: "em_min", ucd: "em.wl;stat.min", utype: "Char.SpectralAxis.Coverage.Bounds.Limits.LoLimit", units: "m" },
    em_max: { name: "em_max", ucd: "em.wl;stat.max", utype: "Char.SpectralAxis.Coverage.Bounds.Limits.HiLimit", units: "m" },
    em_res_power: { name: "em_res_power", ucd: "spect.resolution", utype: "Char.SpectralAxis.Resolution.ResolPower.refVal", units: null },
    em_xel: { name: "em_xel", ucd: "meta.number", utype: "Char.SpectralAxis.numBins", units: null },
    o_ucd: { name: "o_ucd", ucd: "meta.ucd", utype: "Char.ObservableAxis.ucd", units: null },
    pol_states: { name: "pol_states", ucd: "meta.code;phys.polarization", utype: "Char.PolarizationAxis.stateList", units: null },
    pol_xel: { name: "pol_xel", ucd: "meta.number", utype: "Char.PolarizationAxis.numBins", units: null },
    facility_name: { name: "facility_name", ucd: "meta.id;instr.tel", utype: "Provenance.ObsConfig.Facility.name", units: null },
    instrument_name: { name: "instrument_name", ucd: "meta.id;instr", utype: "Provenance.ObsConfig.Instrument.name", units: null }
  }, B.clickOnAccessUrlAction = function(A) {
    hB.parse(A, (g, C) => {
      console.log(g);
    });
  }, B.COLOR = "#004500";
  function B() {
  }
  return B.parseFields = function(A) {
    let g = {};
    const C = B.MANDATORY_FIELDS.s_ra, I = B.MANDATORY_FIELDS.s_dec, E = B.MANDATORY_FIELDS.s_region, o = B.MANDATORY_FIELDS.access_url, t = B.MANDATORY_FIELDS.access_format;
    let s = B.findMandatoryField(A, C.name, C.ucd, C.utype), M = B.findMandatoryField(A, I.name, I.ucd, I.utype), N = B.findMandatoryField(A, E.name, E.ucd, E.utype), k = B.findMandatoryField(A, o.name, o.ucd, o.utype), U = B.findMandatoryField(A, t.name, t.ucd, t.utype), Y = 0;
    return A.forEach((K) => {
      let l = K.name ? K.name : K.id, d;
      Y == s ? d = "s_ra" : Y == M ? d = "s_dec" : Y == N ? d = "s_region" : Y == k ? d = "access_url" : Y == U ? d = "access_format" : d = l, g[d] = {
        name: l,
        idx: Y
      }, Y++;
    }), g;
  }, B.findMandatoryField = function(A, g = null, C = null, I = null) {
    if (Z.isInt(g) && g < A.length)
      return g;
    if (g)
      for (var E = 0, o = A.length; E < o; E++) {
        var t = A[E];
        if (t.ID && t.ID === g || t.name && t.name === g)
          return E;
      }
    if (C)
      for (var s = C.replace(".", "_"), E = 0, o = A.length; E < o; E++) {
        var t = A[E];
        if (t.ucd) {
          var M = $.trim(t.ucd.toLowerCase());
          if (M.indexOf(C) == 0 || M.indexOf(s) == 0)
            return E;
        }
      }
    if (I)
      for (var E = 0, o = A.length; E < o; E++) {
        var t = A[E];
        if (t.utype) {
          var N = $.trim(t.utype.toLowerCase());
          if (N === I)
            return E;
        }
      }
    throw "Mandatory field " + g + " not found";
  }, B.handleActions = function(A) {
    let g = A.fields, C = A.view.aladin;
    A.addFieldClickCallback("access_url", (I) => {
      let E = g.access_url.name, o = g.access_format.name, t = I[E], s = I[o], M = () => {
        var N = g.obs_id && I[g.obs_id.name], k = N || t, U = (K, l, d, q) => {
          C.gotoRaDec(K, l), C.setFoV(d);
        };
        let Y = C.createImageFITS(t, k, {}, U);
        C.setOverlayImageLayer(Y, Z.uuidv4());
      };
      switch (s) {
        case "application/x-votable+xml;content=datalink":
          kE.handleActions(t, C);
          break;
        case "image/fits":
          M();
          break;
        case "application/fits":
          M();
          break;
        case "application/x-fits-mef":
          M();
          break;
        default:
          console.warn("Access format ", s, " not yet implemented or not recognized. Download the file triggered"), Z.download(t);
          break;
      }
    });
  }, B;
}(), hB = function() {
  function B(A, g, C) {
    fetch(A).then((I) => I.text()).then((I) => {
      CA.AL_USE_WASM.dispatchedTo(document.body, { callback: (E) => {
        let o = E.parseVOTable(I);
        g(o);
      } });
    }).catch((I) => C(I));
  }
  return B.parse = function(A, g, C, I, E, o) {
    I && (A = Z.handleCORSNotSameOrigin(A)), fetch(A).then((t) => t.text()).then((t) => {
      CA.AL_USE_WASM.dispatchedTo(document.body, {
        callback: (s) => {
          s.parseVOTable(t).votable.get("resources").forEach((N) => {
            let k = N.get("tables");
            k && k.forEach((U) => {
              let Y = U.get("elems").filter((l) => (l.elem_type || l.get("elem_type")) === "Field").map((l) => Object.fromEntries(l));
              try {
                Y = rQ.parseFields(Y), Y.subtype = "ObsCore";
              } catch {
                Y = bg.parseFields(Y, E, o);
              }
              let K = U.get("data");
              if (K) {
                let l = K.get("rows");
                l && g(Y, l);
              }
            });
          });
        }
      });
    }).catch((t) => {
      if (C)
        C(t);
      else
        throw t;
    });
  }, B;
}(), bg = function() {
  function B(g) {
    g = g || {}, this.uuid = Z.uuidv4(), this.type = "catalog", this.name = g.name || "catalog", this.color = g.color || _A.getNextColor(), this.sourceSize = g.sourceSize || 8, this.markerSize = g.sourceSize || 12, this.selectSize = this.sourceSize, this.shape = g.shape || "square", this.maxNbSources = g.limit || void 0, this.onClick = g.onClick || void 0, this.raField = g.raField || void 0, this.decField = g.decField || void 0, this.filterFn = g.filter || void 0, this.fieldsClickedActions = {}, this.fields = void 0, this.indexationNorder = 5, this.sources = [], this.ra = [], this.dec = [], this.footprints = [], this.displayLabel = g.displayLabel || !1, this.labelColor = g.labelColor || this.color, this.labelFont = g.labelFont || "10px sans-serif", this.displayLabel && (this.labelColumn = g.labelColumn, this.labelColumn || (this.displayLabel = !1)), this.selectionColor = "#00ff00", this.updateShape(g), this.cacheMarkerCanvas = document.createElement("canvas"), this.cacheMarkerCanvas.width = this.markerSize, this.cacheMarkerCanvas.height = this.markerSize;
    var C = this.cacheMarkerCanvas.getContext("2d");
    C.fillStyle = this.color, C.beginPath();
    var I = this.markerSize / 2;
    C.arc(I, I, I - 2, 0, 2 * Math.PI, !1), C.fill(), C.lineWidth = 2, C.strokeStyle = "#ccc", C.stroke(), this.isShowing = !0;
  }
  B.createShape = function(g, C, I) {
    if (g instanceof Image || g instanceof HTMLCanvasElement)
      return g;
    var E = document.createElement("canvas");
    E.width = E.height = I;
    var o = E.getContext("2d");
    return o.beginPath(), o.strokeStyle = C, o.lineWidth = 2, g == "plus" ? (o.moveTo(I / 2, 0), o.lineTo(I / 2, I), o.stroke(), o.moveTo(0, I / 2), o.lineTo(I, I / 2), o.stroke()) : g == "cross" ? (o.moveTo(0, 0), o.lineTo(I - 1, I - 1), o.stroke(), o.moveTo(I - 1, 0), o.lineTo(0, I - 1), o.stroke()) : g == "rhomb" ? (o.moveTo(I / 2, 0), o.lineTo(0, I / 2), o.lineTo(I / 2, I), o.lineTo(I, I / 2), o.lineTo(I / 2, 0), o.stroke()) : g == "triangle" ? (o.moveTo(I / 2, 0), o.lineTo(0, I - 1), o.lineTo(I - 1, I - 1), o.lineTo(I / 2, 0), o.stroke()) : g == "circle" ? (o.arc(I / 2, I / 2, I / 2 - 1, 0, 2 * Math.PI, !0), o.stroke()) : (o.moveTo(1, 0), o.lineTo(1, I - 1), o.lineTo(I - 1, I - 1), o.lineTo(I - 1, 1), o.lineTo(1, 1), o.stroke()), E;
  };
  function A(g, C, I) {
    var E, o;
    if (E = o = null, C)
      for (var t = 0, s = g.length; t < s; t++) {
        var M = g[t];
        if (Z.isInt(C) && C < g.length) {
          E = C;
          break;
        }
        if (M.ID && M.ID === C || M.name && M.name === C) {
          E = t;
          break;
        }
      }
    if (I)
      for (var t = 0, s = g.length; t < s; t++) {
        var M = g[t];
        if (Z.isInt(I) && I < g.length) {
          o = I;
          break;
        }
        if (M.ID && M.ID === I || M.name && M.name === I) {
          o = t;
          break;
        }
      }
    for (var t = 0, s = g.length; t < s && !(E != null && o != null); t++) {
      var M = g[t];
      if (!E && M.ucd) {
        var N = O.trim(M.ucd.toLowerCase());
        if (N.indexOf("pos.eq.ra") == 0 || N.indexOf("pos_eq_ra") == 0) {
          E = t;
          continue;
        }
      }
      if (!o && M.ucd) {
        var N = O.trim(M.ucd.toLowerCase());
        if (N.indexOf("pos.eq.dec") == 0 || N.indexOf("pos_eq_dec") == 0) {
          o = t;
          continue;
        }
      }
    }
    if (E == null && o == null)
      for (var t = 0, s = g.length; t < s; t++) {
        var M = g[t], k = M.name || M.ID || "";
        if (k = k.toLowerCase(), !E && (k.indexOf("ra") == 0 || k.indexOf("_ra") == 0 || k.indexOf("ra(icrs)") == 0 || k.indexOf("_ra") == 0 || k.indexOf("alpha") == 0)) {
          E = t;
          continue;
        }
        if (!o && (k.indexOf("dej2000") == 0 || k.indexOf("_dej2000") == 0 || k.indexOf("de") == 0 || k.indexOf("de(icrs)") == 0 || k.indexOf("_de") == 0 || k.indexOf("delta") == 0)) {
          o = t;
          continue;
        }
      }
    return (E == null || o == null) && (E = 0, o = 1), [E, o];
  }
  return B.parseFields = function(g, C, I) {
    let [E, o] = A(g, C, I), t = {}, s = 0;
    return g.forEach((M) => {
      let N = M.name ? M.name : M.id, k;
      s == E ? k = "ra" : s == o ? k = "dec" : k = N, t[k] = {
        name: N,
        idx: s
      }, s++;
    }), t;
  }, B.parseVOTable = function(g, C, I, E, o, t, s) {
    hB.parse(
      g,
      (M, N) => {
        let k = [], U = [];
        var Y = new $A();
        N.every((K) => {
          let l, d, q;
          var W = {};
          for (const [a, iA] of Object.entries(M)) {
            a === "s_region" ? q = K[iA.idx] : a === "ra" || a === "s_ra" ? l = K[iA.idx] : (a === "dec" || a === "s_dec") && (d = K[iA.idx]);
            var T = iA.name;
            W[T] = K[iA.idx];
          }
          let v = null;
          l && d && ((!Z.isNumber(l) || !Z.isNumber(d)) && (Y.parse(l + " " + d), l = Y.lon, d = Y.lat), v = new GB(l, d, W));
          let _ = null;
          if (q) {
            let a = zA.footprintsFromSTCS(q, { lineWidth: 2 });
            _ = new yQ(a, v);
          }
          if (_)
            U.push(_);
          else if (v && (k.push(v), E && k.length == E))
            return !1;
          return !0;
        }), C && C(k, U, M);
      },
      I,
      o,
      t,
      s
    );
  }, B.prototype.updateShape = function(g) {
    g = g || {}, this.color = g.color || this.color || _A.getNextColor(), this.sourceSize = g.sourceSize || this.sourceSize || 6, this.shape = g.shape || this.shape || "square", this._shapeIsFunction = !1, typeof this.shape == "function" && (this._shapeIsFunction = !0), (this.shape instanceof Image || this.shape instanceof HTMLCanvasElement) && (this.sourceSize = this.shape.width), this.selectSize = this.sourceSize + 2, this.cacheCanvas = B.createShape(this.shape, this.color, this.sourceSize), this.cacheSelectCanvas = B.createShape(this.shape, this.selectionColor, this.selectSize), this.cacheHoverCanvas = B.createShape(this.shape, this.hoverColor, this.sourceSize), this.reportChange();
  }, B.prototype.addSources = function(g) {
    if (g = [].concat(g), g.length !== 0) {
      if (!this.fields) {
        let o = [];
        for (var C in g[0].data)
          o.push({ name: C });
        o = B.parseFields(o, this.raField, this.decField), this.setFields(o);
      }
      this.sources = this.sources.concat(g);
      for (var I = 0, E = g.length; I < E; I++)
        g[I].setCatalog(this), this.ra.push(g[I].ra), this.dec.push(g[I].dec);
      CA.AL_USE_WASM.dispatchedTo(document.body, {
        callback: (o) => {
          o.setCatalog(this), this.reportChange();
        }
      });
    }
  }, B.prototype.addFootprints = function(g) {
    g = [].concat(g), this.footprints = this.footprints.concat(g);
    for (var C = 0, I = g.length; C < I; C++)
      g[C].setCatalog(this), g[C].setColor(this.color), g[C].setSelectionColor(this.selectionColor);
    this.reportChange();
  }, B.prototype.setFields = function(g) {
    this.fields = g;
  }, B.prototype.addFieldClickCallback = function(g, C) {
    this.fieldsClickedActions[g] = C;
  }, B.prototype.isObsCore = function() {
    return this.fields && this.fields.subtype === "ObsCore";
  }, B.prototype.addSourcesAsArray = function(g, C) {
    for (var I = [], E = 0; E < g.length; E++)
      I.push({ name: g[E] });
    I = B.parseFields(I, this.raField, this.decField), this.setFields(I);
    var o, t;
    o = I.ra.idx, t = I.dec.idx;
    for (var s = [], M = new $A(), N, k, U, Y, K = 0; K < C.length; K++) {
      U = C[K], Z.isNumber(U[o]) && Z.isNumber(U[t]) ? (N = parseFloat(U[o]), k = parseFloat(U[t])) : (M.parse(U[o] + " " + U[t]), N = M.lon, k = M.lat), Y = {};
      for (var E = 0; E < g.length; E++)
        Y[g[E]] = U[E];
      s.push(zA.source(N, k, Y));
    }
    this.addSources(s);
  }, B.prototype.getSources = function() {
    return this.sources;
  }, B.prototype.getFootprints = function() {
    return this.footprints;
  }, B.prototype.selectAll = function() {
    if (this.sources)
      for (var g = 0; g < this.sources.length; g++)
        this.sources[g].select();
  }, B.prototype.deselectAll = function() {
    if (this.sources)
      for (var g = 0; g < this.sources.length; g++)
        this.sources[g].deselect();
  }, B.prototype.getSource = function(g) {
    return g < this.sources.length ? this.sources[g] : null;
  }, B.prototype.setView = function(g) {
    this.view = g, this.reportChange();
  }, B.prototype.setColor = function(g) {
    this.color = g, this.updateShape();
  }, B.prototype.setSelectionColor = function(g) {
    this.selectionColor = g, this.updateShape();
  }, B.prototype.setHoverColor = function(g) {
    this.hoverColor = g, this.updateShape();
  }, B.prototype.setSourceSize = function(g) {
    this.sourceSize = g, this.updateShape();
  }, B.prototype.setShape = function(g) {
    this.shape = g, this.updateShape();
  }, B.prototype.getSourceSize = function() {
    return this.sourceSize;
  }, B.prototype.remove = function(g) {
    var C = this.sources.indexOf(g);
    C < 0 || (this.sources[C].deselect(), this.sources.splice(C, 1), this.ra.splice(C, 1), this.dec.splice(C, 1), this.reportChange());
  }, B.prototype.removeAll = B.prototype.clear = function() {
    this.sources = [], this.ra = [], this.dec = [], this.footprints = [];
  }, B.prototype.draw = function(g, C, I, E, o, t) {
    if (!this.isShowing)
      return;
    g.strokeStyle = this.color, this._shapeIsFunction && g.save();
    const s = this.drawSources(g, I, E);
    this._shapeIsFunction && g.restore(), this.displayLabel && (g.fillStyle = this.labelColor, g.font = this.labelFont, s.forEach((M) => {
      this.drawSourceLabel(M, g);
    })), this.drawFootprints(g);
  }, B.prototype.drawSources = function(g, C, I) {
    if (!this.sources)
      return;
    let E = [], o = this.view.wasm.worldToScreenVec(this.ra, this.dec), t = this;
    return this.sources.forEach(function(s, M) {
      o[2 * M] && o[2 * M + 1] && (!t.filterFn || t.filterFn(s)) && (s.x = o[2 * M], s.y = o[2 * M + 1], t.drawSource(s, g, C, I), E.push(s));
    }), E;
  }, B.prototype.drawSource = function(g, C, I, E) {
    return g.isShowing && g.x <= I && g.x >= 0 && g.y <= E && g.y >= 0 ? (this._shapeIsFunction ? this.shape(g, C, this.view.getViewParams()) : g.marker && g.useMarkerDefaultIcon ? C.drawImage(this.cacheMarkerCanvas, g.x - this.sourceSize / 2, g.y - this.sourceSize / 2) : g.isSelected ? C.drawImage(this.cacheSelectCanvas, g.x - this.selectSize / 2, g.y - this.selectSize / 2) : C.drawImage(this.cacheCanvas, g.x - this.cacheCanvas.width / 2, g.y - this.cacheCanvas.height / 2), g.popup && g.popup.setPosition(g.x, g.y), !0) : !1;
  }, B.prototype.drawSourceLabel = function(g, C) {
    if (!(!g || !g.isShowing || !g.x || !g.y)) {
      var I = g.data[this.labelColumn];
      I && C.fillText(I, g.x, g.y);
    }
  }, B.prototype.drawFootprints = function(g) {
    this.footprints.forEach((C) => {
      C.draw(g, this.view);
    });
  }, B.prototype.reportChange = function() {
    this.view && this.view.requestRedraw();
  }, B.prototype.show = function() {
    this.isShowing || (this.isShowing = !0, this.footprints && this.footprints.forEach((g) => g.show()), this.reportChange());
  }, B.prototype.hide = function() {
    this.isShowing && (this.isShowing = !1, this.view && this.view.popup && this.view.popup.source && this.view.popup.source.catalog == this && this.view.popup.hide(), this.footprints && this.footprints.forEach((g) => g.hide()), this.reportChange());
  }, B;
}(), NQ = function() {
  let B = function(C, I, E, o) {
    o = o || {}, this.uuid = Z.uuidv4(), this.type = "progressivecat", this.rootUrl = C, Z.isHttpsContext() && (/u-strasbg.fr/i.test(this.rootUrl) || /unistra.fr/i.test(this.rootUrl)) && (this.rootUrl = this.rootUrl.replace("http://", "https://")), this.frameStr = I, this.frame = RA.fromString(I) || RA.J2000, this.maxOrder = parseInt(E), this.isShowing = !0, this.name = o.name || "progressive-cat", this.color = o.color || _A.getNextColor(), this.shape = o.shape || "square", this.sourceSize = o.sourceSize || 6, this.selectSize = this.sourceSize + 2, this.selectionColor = "#00ff00", this.filterFn = o.filter || void 0, this.onClick = o.onClick || void 0, this.sourcesCache = new Z.LRUCache(256), (this.shape instanceof Image || this.shape instanceof HTMLCanvasElement) && (this.sourceSize = this.shape.width), this._shapeIsFunction = !1, typeof this.shape == "function" && (this._shapeIsFunction = !0), this.updateShape(o), this.maxOrderAllsky = 2, this.isReady = !1, this.tilesInView = [];
  };
  B.readProperties = function(C, I, E) {
    if (I) {
      var o = C + "/properties";
      O.ajax({
        url: o,
        method: "GET",
        dataType: "text",
        success: function(t) {
          for (var s = {}, M = t.split(`
`), N = 0; N < M.length; N++) {
            var k = M[N], U = k.indexOf("="), Y = O.trim(k.substring(0, U)), K = O.trim(k.substring(U + 1));
            s[Y] = K;
          }
          I(s);
        },
        error: function(t) {
          E && E(t);
        }
      });
    }
  };
  function A(C, I) {
    var E = ["name", "ID", "ucd", "utype", "unit", "datatype", "arraysize", "width", "precision"], o = [], t = 0;
    return C.keyRa = C.keyDec = null, O(I).find("FIELD").each(function() {
      for (var s = {}, M = 0; M < E.length; M++) {
        var N = E[M];
        O(this).attr(N) && (s[N] = O(this).attr(N));
      }
      s.ID || (s.ID = "col_" + t), !C.keyRa && s.ucd && (s.ucd.indexOf("pos.eq.ra") == 0 || s.ucd.indexOf("POS_EQ_RA") == 0) && (s.name ? C.keyRa = s.name : C.keyRa = s.ID), !C.keyDec && s.ucd && (s.ucd.indexOf("pos.eq.dec") == 0 || s.ucd.indexOf("POS_EQ_DEC") == 0) && (s.name ? C.keyDec = s.name : C.keyDec = s.ID), o.push(s), t++;
    }), o;
  }
  function g(C, I, E) {
    if (!C.keyRa || !C.keyDec)
      return [];
    for (var o = I.split(`
`), t = [], s = 0; s < E.length; s++)
      E[s].name ? t.push(E[s].name) : t.push(E[s].ID);
    for (var M = [], N = new $A(), k, U = 2; U < o.length; U++) {
      var Y = {}, K = o[U].split("	");
      if (!(K.length < t.length)) {
        for (var l = 0; l < t.length; l++)
          Y[t[l]] = K[l];
        var d, q;
        Z.isNumber(Y[C.keyRa]) && Z.isNumber(Y[C.keyDec]) ? (d = parseFloat(Y[C.keyRa]), q = parseFloat(Y[C.keyDec])) : (N.parse(Y[C.keyRa] + " " + Y[C.keyDec]), d = N.lon, q = N.lat), k = new GB(d, q, Y), M.push(k), k.setCatalog(C);
      }
    }
    return M;
  }
  return B.prototype = {
    init: function(C) {
      var I = this;
      this.view = C, this.maxOrder && this.frameStr ? this._loadMetadata() : B.readProperties(
        I.rootUrl,
        function(E) {
          I.properties = E, I.maxOrder = parseInt(I.properties.hips_order), I.frame = RA.fromString(I.properties.hips_frame), I._loadMetadata();
        },
        function(E) {
          console.log("Could not find properties for HiPS " + I.rootUrl);
        }
      );
    },
    updateShape: bg.prototype.updateShape,
    _loadMetadata: function() {
      var C = this;
      O.ajax({
        url: C.rootUrl + "/Metadata.xml",
        method: "GET",
        success: function(I) {
          C.fields = A(C, I), C._loadAllskyNewMethod();
        },
        error: function(I) {
          C._loadAllskyOldMethod();
        }
      });
    },
    _loadAllskyNewMethod: function() {
      var C = this;
      O.ajax({
        url: C.rootUrl + "/Norder1/Allsky.tsv",
        method: "GET",
        success: function(I) {
          C.order1Sources = g(C, I, C.fields), C.order2Sources && (C.isReady = !0, C._finishInitWhenReady());
        },
        error: function(I) {
          console.log("Something went wrong: " + I);
        }
      }), O.ajax({
        url: C.rootUrl + "/Norder2/Allsky.tsv",
        method: "GET",
        success: function(I) {
          C.order2Sources = g(C, I, C.fields), C.order1Sources && (C.isReady = !0, C._finishInitWhenReady());
        },
        error: function(I) {
          console.log("Something went wrong: " + I);
        }
      });
    },
    _loadAllskyOldMethod: function() {
      this.maxOrderAllsky = 3, this._loadLevel2Sources(), this._loadLevel3Sources();
    },
    _loadLevel2Sources: function() {
      var C = this;
      O.ajax({
        url: C.rootUrl + "/Norder2/Allsky.xml",
        method: "GET",
        success: function(I) {
          C.fields = A(C, I), C.order2Sources = g(C, O(I).find("CSV").text(), C.fields), C.order3Sources && (C.isReady = !0, C._finishInitWhenReady());
        },
        error: function(I) {
          console.log("Something went wrong: " + I);
        }
      });
    },
    _loadLevel3Sources: function() {
      var C = this;
      O.ajax({
        url: C.rootUrl + "/Norder3/Allsky.xml",
        method: "GET",
        success: function(I) {
          C.order3Sources = g(C, O(I).find("CSV").text(), C.fields), C.order2Sources && (C.isReady = !0, C._finishInitWhenReady());
        },
        error: function(I) {
          console.log("Something went wrong: " + I);
        }
      });
    },
    _finishInitWhenReady: function() {
      this.view.requestRedraw(), this.loadNeededTiles();
    },
    draw: function(C, I, E, o, t, s) {
      if (!this.isShowing || !this.isReady)
        return;
      this._shapeIsFunction && C.save(), this.drawSources(this.order1Sources, C, E, o), this.view.realNorder >= 1 && this.drawSources(this.order2Sources, C, E, o), this.maxOrderAllsky === 3 && this.view.realNorder >= 2 && this.drawSources(this.order3Sources, C, E, o);
      let M, N;
      this.tilesInView.forEach((k) => {
        M = k[0] + "-" + k[1], N = this.sourcesCache.get(M), N && this.drawSources(N, C, E, o);
      }), this._shapeIsFunction && C.restore();
    },
    drawSources: function(C, I, E, o) {
      if (!C)
        return;
      let t = [], s = [];
      C.forEach((k) => {
        t.push(k.ra), s.push(k.dec);
      });
      let M = this.view.wasm.worldToScreenVec(t, s), N = this;
      C.forEach(function(k, U) {
        M[2 * U] && M[2 * U + 1] && (!N.filterFn || N.filterFn(k)) && (k.x = M[2 * U], k.y = M[2 * U + 1], N.drawSource(k, I, E, o));
      });
    },
    drawSource: bg.prototype.drawSource,
    getSources: function() {
      var C = [];
      if (this.order1Sources && (C = C.concat(this.order1Sources)), this.order2Sources && (C = C.concat(this.order2Sources)), this.order3Sources && (C = C.concat(this.order3Sources)), this.tilesInView)
        for (var I, E, o, t = 0; t < this.tilesInView.length; t++)
          o = this.tilesInView[t], E = o[0] + "-" + o[1], I = this.sourcesCache.get(E), I && (C = C.concat(I));
      return C;
    },
    deselectAll: function() {
      if (this.order1Sources)
        for (var C = 0; C < this.order1Sources.length; C++)
          this.order1Sources[C].deselect();
      if (this.order2Sources)
        for (var C = 0; C < this.order2Sources.length; C++)
          this.order2Sources[C].deselect();
      if (this.order3Sources)
        for (var C = 0; C < this.order3Sources.length; C++)
          this.order3Sources[C].deselect();
      var I = this.sourcesCache.keys();
      for (key in I)
        if (this.sourcesCache[key])
          for (var E = this.sourcesCache[key], C = 0; C < E.length; C++)
            E[C].deselect();
    },
    show: function() {
      this.isShowing || (this.isShowing = !0, this.loadNeededTiles(), this.reportChange());
    },
    hide: function() {
      this.isShowing && (this.isShowing = !1, this.reportChange());
    },
    reportChange: function() {
      this.view.requestRedraw();
    },
    getTileURL: function(C, I) {
      var E = Math.floor(I / 1e4) * 1e4;
      return this.rootUrl + "/Norder" + C + "/Dir" + E + "/Npix" + I + ".tsv";
    },
    // todo, allow HiPS cats to support footprints 
    getFootprints: function() {
      return null;
    },
    loadNeededTiles: function() {
      if (!this.isShowing)
        return;
      this.tilesInView = [];
      var C = this.view.realNorder;
      C > this.maxOrder && (C = this.maxOrder);
      var I = this.view.getVisibleCells(C);
      let E = C;
      for (; I.length > 12 && E > this.maxOrderAllsky; )
        E--, I = this.view.getVisibleCells(E);
      if (C = E, !(C <= this.maxOrderAllsky)) {
        for (var o, t, s = this.maxOrderAllsky + 1; s <= C; s++) {
          o = [];
          for (var M = 0; M < I.length; M++)
            t = Math.floor(I[M].ipix / Math.pow(4, C - s)), o.indexOf(t) < 0 && o.push(t);
          for (var N = 0; N < o.length; N++)
            this.tilesInView.push([s, o[N]]);
        }
        for (var k, U, M = 0; M < this.tilesInView.length; M++)
          k = this.tilesInView[M], U = k[0] + "-" + k[1], this.sourcesCache.get(U) || function(K, l, d) {
            var q = l + "-" + d;
            O.ajax({
              /*
              url: Aladin.JSONP_PROXY,
              data: {"url": self.getTileURL(norder, ipix)},
              */
              // ATTENTIOn : je passe en JSON direct, car je n'arrive pas a choper les 404 en JSONP
              url: K.getTileURL(l, d),
              method: "GET",
              //dataType: 'jsonp',
              success: function(W) {
                K.sourcesCache.set(q, g(K, W, K.fields)), K.view.requestRedraw();
              },
              error: function() {
                K.sourcesCache.set(q, []);
              }
            });
          }(this, k[0], k[1]);
      }
    },
    reportChange: function() {
      this.view && this.view.requestRedraw();
    }
  }, B;
}(), nE = function() {
  let B = {};
  return B.cache = {}, B.SESAME_URL = "http://cds.u-strasbg.fr/cgi-bin/nph-sesame.jsonp", B.getTargetRADec = function(A, g, C) {
    if (g) {
      var I = /[a-zA-Z]/.test(A);
      if (I)
        B.resolve(
          A,
          function(o) {
            g({
              ra: o.Target.Resolver.jradeg,
              dec: o.Target.Resolver.jdedeg
            });
          },
          function(o) {
            C && C();
          }
        );
      else {
        var E = new Coo();
        E.parse(A), g && g({ ra: E.lon, dec: E.lat });
      }
    }
  }, B.resolve = function(A, g, C) {
    var I = B.SESAME_URL;
    Z.isHttpsContext() && (I = I.replace("http://", "https://")), O.ajax({
      url: I,
      data: { object: A },
      method: "GET",
      dataType: "jsonp",
      success: function(E) {
        E.Target && E.Target.Resolver && E.Target.Resolver ? g(E) : C(E);
      },
      error: C
    });
  }, B;
}(), FE = function() {
  let B = {};
  B.cache = {}, B.URL = "https://alasky.cds.unistra.fr/planetary-features/resolve";
  function A(g) {
    let C = "", I = [""], E = [I], o = 0, t = 0, s = !0, M;
    for (M of g)
      M === '"' ? (s && M === C && (I[o] += M), s = !s) : M === "," && s ? M = I[++o] = "" : M === `
` && s ? (C === "\r" && (I[o] = I[o].slice(0, -1)), I = E[++t] = [M = ""], o = 0) : I[o] += M, C = M;
    return E;
  }
  return B.resolve = function(g, C, I, E) {
    const o = B.URL;
    O.ajax({
      url: o,
      data: { identifier: g, body: C, threshold: 0.7, format: "csv" },
      method: "GET",
      dataType: "text",
      success: function(t) {
        const s = t.split(`
`), M = A(s[0])[0];
        if (s.length > 1 && s[1].length > 0) {
          const N = A(s[1])[0], k = M.findIndex((q) => q.includes("longitude")), U = M.findIndex((q) => q.includes("latitude"));
          let Y = parseFloat(N[k]);
          const K = parseFloat(N[U]);
          let l = !0;
          const d = M.indexOf("coordinate_system");
          d > 0 && N[d].includes("+West") && (l = !1), l || (Y = 360 - Y), I({ lon: Y, lat: K });
        } else
          E(data);
      },
      error: E
    });
  }, B;
}(), RE = function() {
  function B(A) {
    this.isShowing = !1;
    let g = document.createElement("div");
    g.setAttribute("class", "aladin-measurement-div"), this.element = g, this.savedTablesIdx = 0, this.savedTables = [], A.appendChild(this.element);
  }
  return B.prototype.updateRows = function() {
    let A = this.element.querySelector("tbody");
    A.innerHTML = "";
    let g = this.tables[this.curTableIdx], C = "";
    if (g.rows.forEach((I) => {
      C += "<tr>";
      for (let E in I.data) {
        const o = I.data[E] || "--";
        if (C += '<td class="' + E + '">', typeof o == "string")
          try {
            let t = new URL(o), s = "<a href=" + t + ' target="_blank">' + t + "</a>";
            C += s;
          } catch {
            C += o;
          }
        else
          C += o;
        C += "</td>";
      }
      C += "</tr>";
    }), A.innerHTML = C, g.fieldsClickedActions)
      for (let I in g.fieldsClickedActions)
        A.querySelectorAll("." + I).forEach(function(E, o) {
          E.addEventListener("click", (t) => {
            let s = g.fieldsClickedActions[I];
            s(g.rows[o].data), t.preventDefault();
          }, !1);
        });
  }, B.prototype.showPreviousMeasurement = function() {
    this.savedTablesIdx--, this.savedTablesIdx < 0 && (this.savedTablesIdx = 0);
    let A = this.savedTables[this.savedTablesIdx];
    A && (this.update(A), this.updateStateNavigation());
  }, B.prototype.showNextMeasurement = function() {
    this.savedTablesIdx++, this.savedTablesIdx >= this.savedTables.length && (this.savedTablesIdx = this.savedTables.length - 1);
    let A = this.savedTables[this.savedTablesIdx];
    A && (this.update(A), this.updateStateNavigation());
  }, B.prototype.showMeasurement = function(A, g) {
    A.length !== 0 && (this.update(A), g && g.save && (this.saveState(), this.updateStateNavigation()));
  }, B.prototype.updateStateNavigation = function() {
    let A = this.element.querySelector(".tabs");
    if (this.savedTables.length >= 2) {
      let g = document.createElement("button");
      g.setAttribute("title", "Go back to the previous table"), this.savedTablesIdx == 0 && (g.disabled = !0), g.addEventListener(
        "click",
        () => this.showPreviousMeasurement(),
        !1
      ), g.innerText = "<", A.appendChild(g);
      let C = document.createElement("button");
      C.setAttribute("title", "Go to the next table"), (this.savedTables.length == 0 || this.savedTablesIdx == this.savedTables.length - 1) && (C.disabled = !0), C.addEventListener(
        "click",
        () => this.showNextMeasurement(),
        !1
      ), C.innerText = ">", A.appendChild(C);
    }
  }, B.prototype.saveState = function() {
    this.savedTables.length === 0 ? this.savedTables.push(this.tables) : this.tables !== this.savedTables[this.savedTablesIdx] && (this.savedTables = this.savedTables.slice(0, this.savedTablesIdx + 1), this.savedTables.push(this.tables), this.savedTablesIdx = this.savedTables.length - 1);
  }, B.prototype.update = function(A) {
    this.tables = A, this.curTableIdx = 0;
    let g = A[this.curTableIdx];
    this.element.innerHTML = "";
    let C = this.createTabs();
    this.element.appendChild(C);
    let I = document.createElement("table");
    I.style.borderColor = g.color;
    const E = B.createTableHeader(g), o = document.createElement("tbody");
    I.appendChild(E), I.appendChild(o), this.element.appendChild(I), this.updateRows(), this.show();
  }, B.prototype.createTabs = function() {
    let A = document.createElement("div");
    A.setAttribute("class", "tabs");
    let g = this;
    return this.tables.forEach(function(C, I) {
      let E = document.createElement("button");
      E.setAttribute("title", C.name), E.innerText = C.name, E.style.overflow = "hidden", E.style.textOverflow = "ellipsis", E.style.whiteSpace = "nowrap", E.style.maxWidth = "20%", E.addEventListener(
        "click",
        () => {
          g.curTableIdx = I;
          let M = g.element.querySelector("table");
          M.style.borderColor = C.color;
          let N = g.element.querySelector("thead");
          N.parentNode.replaceChild(B.createTableHeader(C), N), g.updateRows();
        },
        !1
      ), E.style.backgroundColor = C.color;
      let o = _A.standardizeColor(C.color), t = _A.hexToRgb(o);
      t = "rgb(" + t.r + ", " + t.g + ", " + t.b + ")";
      let s = _A.getLabelColorForBackground(t);
      E.style.color = s, A.appendChild(E);
    }), A;
  }, B.createTableHeader = function(A) {
    let g = document.createElement("thead");
    var C = "<tr>";
    for (let [I, E] of Object.entries(A.fields))
      C += "<th>" + E.name + "</th>";
    return C += "</thead>", g.innerHTML = C, g;
  }, B.prototype.show = function() {
    this.element.style.visibility = "visible";
  }, B.prototype.hide = function() {
    this.savedTables = [], this.savedTablesIdx = 0, this.curTableIdx = 0, this.element.style.visibility = "hidden";
  }, B;
}(), KE = function() {
  function B(A) {
    this.$div = O(A);
  }
  return B.prototype.update = function(A, g, C, I) {
    var E = new $A(A, g, 7);
    C == RA.J2000 ? this.$div.html(E.format("s/")) : C == RA.J2000d ? this.$div.html(E.format("d/")) : this.$div.html(E.format("d/")), this.$div.toggleClass("aladin-reticleColor", I);
  }, B;
}(), Fg = {};
Fg.update = function(B) {
  const A = Fg.LAYERS.find(({ id: C }) => B.id.endsWith(C)), g = {
    ...B.colorCfg.get(),
    imgFormat: B.imgFormat,
    longitudeReversed: B.longitudeReversed
  };
  A ? A.options = g : Fg.LAYERS.push({
    id: B.id,
    name: B.name,
    url: B.url,
    options: g,
    subtype: B.subtype
  });
};
Fg.LAYERS = [
  {
    id: "P/2MASS/color",
    name: "2MASS colored",
    url: "https://alasky.cds.unistra.fr/2MASS/Color",
    maxOrder: 9,
    subtype: "survey"
  },
  {
    id: "P/DSS2/color",
    name: "DSS colored",
    url: "https://alasky.cds.unistra.fr/DSS/DSSColor",
    maxOrder: 9,
    subtype: "survey"
  },
  {
    id: "P/DSS2/red",
    name: "DSS2 Red (F+R)",
    url: "https://alasky.cds.unistra.fr/DSS/DSS2Merged",
    maxOrder: 9,
    subtype: "survey",
    // options
    options: {
      minCut: 1e3,
      maxCut: 1e4,
      colormap: "magma",
      stretch: "Linear",
      imgFormat: "fits"
    }
  },
  {
    id: "P/DM/I/350/gaiaedr3",
    name: "Density map for Gaia EDR3 (I/350/gaiaedr3)",
    url: "https://alasky.cds.unistra.fr/ancillary/GaiaEDR3/density-map",
    maxOrder: 7,
    subtype: "survey",
    // options
    options: {
      minCut: 0,
      maxCut: 12e3,
      stretch: "asinh",
      colormap: "rdylbu",
      imgFormat: "fits"
    }
  },
  {
    id: "P/PanSTARRS/DR1/g",
    name: "PanSTARRS DR1 g",
    url: "https://alasky.cds.unistra.fr/Pan-STARRS/DR1/g",
    maxOrder: 11,
    subtype: "survey",
    // options
    options: {
      minCut: -34,
      maxCut: 7e3,
      stretch: "asinh",
      colormap: "redtemperature",
      imgFormat: "fits"
    }
  },
  {
    id: "P/PanSTARRS/DR1/color-z-zg-g",
    name: "PanSTARRS DR1 color",
    url: "https://alasky.cds.unistra.fr/Pan-STARRS/DR1/color-z-zg-g",
    maxOrder: 11,
    subtype: "survey"
  },
  {
    id: "P/DECaPS/DR1/color",
    name: "DECaPS DR1 color",
    url: "https://alasky.cds.unistra.fr/DECaPS/DR1/color",
    maxOrder: 11,
    subtype: "survey"
  },
  {
    id: "P/Fermi/color",
    name: "Fermi color",
    url: "https://alasky.cds.unistra.fr/Fermi/Color",
    maxOrder: 3,
    subtype: "survey"
  },
  {
    id: "P/Finkbeiner",
    name: "Halpha",
    url: "https://alasky.cds.unistra.fr/FinkbeinerHalpha",
    maxOrder: 3,
    subtype: "survey",
    // options
    options: {
      minCut: -10,
      maxCut: 800,
      colormap: "rdbu",
      imgFormat: "fits"
    }
  },
  {
    id: "P/GALEXGR6_7/NUV",
    name: "GALEXGR6_7 NUV",
    url: "http://alasky.cds.unistra.fr/GALEX/GALEXGR6_7_NUV/",
    maxOrder: 8,
    subtype: "survey"
  },
  {
    id: "P/IRIS/color",
    name: "IRIS colored",
    url: "https://alasky.cds.unistra.fr/IRISColor",
    maxOrder: 3,
    subtype: "survey"
  },
  {
    id: "P/Mellinger/color",
    name: "Mellinger colored",
    url: "https://alasky.cds.unistra.fr/MellingerRGB",
    maxOrder: 4,
    subtype: "survey"
  },
  {
    id: "P/SDSS9/color",
    name: "SDSS9 colored",
    url: "https://alasky.cds.unistra.fr/SDSS/DR9/color",
    maxOrder: 10,
    subtype: "survey"
  },
  {
    id: "P/SDSS9/g",
    name: "SDSS9 band-g",
    url: "https://alasky.cds.unistra.fr/SDSS/DR9/band-g",
    maxOrder: 10,
    subtype: "survey",
    options: {
      stretch: "asinh",
      colormap: "redtemperature",
      imgFormat: "fits"
    }
  },
  {
    id: "P/SPITZER/color",
    name: "IRAC color I1,I2,I4 - (GLIMPSE, SAGE, SAGE-SMC, SINGS)",
    url: "http://alasky.cds.unistra.fr/Spitzer/SpitzerI1I2I4color/",
    maxOrder: 9,
    subtype: "survey"
  },
  {
    id: "P/VTSS/Ha",
    name: "VTSS-Ha",
    url: "https://alasky.cds.unistra.fr/VTSS/Ha",
    maxOrder: 3,
    subtype: "survey",
    options: {
      minCut: -10,
      maxCut: 100,
      colormap: "grayscale",
      imgFormat: "fits"
    }
  },
  {
    id: "xcatdb/P/XMM/PN/color",
    name: "XMM PN colored",
    url: "https://alasky.cds.unistra.fr/cgi/JSONProxy?url=https://saada.unistra.fr/PNColor",
    maxOrder: 7,
    subtype: "survey"
  },
  {
    id: "P/allWISE/color",
    name: "AllWISE color",
    url: "https://alasky.cds.unistra.fr/AllWISE/RGB-W4-W2-W1/",
    maxOrder: 8,
    subtype: "survey"
  },
  {
    id: "P/GLIMPSE360",
    name: "GLIMPSE360",
    // This domain is not giving the CORS headers
    // We need to query by with a proxy equipped with CORS header.
    url: "https://alasky.cds.unistra.fr/cgi/JSONProxy?url=https://www.spitzer.caltech.edu/glimpse360/aladin/data",
    subtype: "survey",
    options: {
      maxOrder: 9,
      imgFormat: "jpg",
      minOrder: 3
    }
  }
];
Fg.getAvailableLayers = function() {
  return Fg.LAYERS;
};
let kQ = function() {
  function B(A) {
    this.properties = A, this.id = this.getID(), this.obsTitle = A.obs_title, this.frame = A.hips_frame, this.order = parseInt(A.hips_order), this.clientSortKey = A.client_sort_key, this.tileFormats = A.hasOwnProperty("hips_tile_format") && A.hips_tile_format.split(" "), this.urls = [], this.urls.push(A.hips_service_url);
    for (var g = 1; A.hasOwnProperty("hips_service_url_" + g); )
      this.urls.push(A["hips_service_url_" + g]), g++;
    this.clientApplications = A.client_application;
  }
  return B.prototype = {
    getServiceURLs: function(A) {
    },
    // return the ID according to the properties
    getID: function() {
      if (this.properties.hasOwnProperty("ID"))
        return this.properties.ID;
      var A = null;
      return this.properties.hasOwnProperty("creator_did") && (A = this.properties.creator_did), A == null && this.properties.hasOwnProperty("publisher_did") && (A = this.properties.publisher_did), A != null && (A.slice(0, 6) === "ivo://" && (A = A.slice(6)), A = A.replace(/\?/g, "/")), A;
    }
  }, B.parseHiPSProperties = function(A) {
    if (A == null)
      return null;
    var g = {};
    A = A.replace(/[\r]/g, "");
    for (var C = A.split(`
`), I = 0; I < C.length; I++) {
      var E = O.trim(C[I]);
      if (E.slice(0, 1) !== "#") {
        var o = E.indexOf("=");
        if (!(o < 0)) {
          var t = O.trim(E.slice(0, o)), s = O.trim(E.slice(o + 1));
          g[t] = s;
        }
      }
    }
    return g;
  }, B.fromURL = function(A, g) {
    var C, I;
    A.slice(-10) === "properties" ? (I = A, C = I.slice(0, -11)) : (A.slice(-1) === "/" && (A = A.slice(0, -1)), C = A, I = C + "/properties");
    var E = function(t) {
      var s = B.parseHiPSProperties(t);
      s.hasOwnProperty("hips_service_url") || (s.hips_service_url = C), typeof g == "function" && g(new B(s));
    }, o = Z.getAjaxObject(I, "GET", "text", !1);
    o.done(function(t) {
      E(t);
    }).fail(function() {
      var t = Z.getAjaxObject(I, "GET", "text", !0);
      t.done(function(s) {
        E(s);
      }).fail(function() {
        typeof g == "function" && g(null);
      });
    });
  }, B.fromProperties = function(A) {
    return new B(A);
  }, B;
}();
class Lg {
  static MIRRORS_HTTP = [
    "http://alaskybis.unistra.fr/MocServer/query",
    "http://alasky.unistra.fr/MocServer/query"
  ];
  // list of base URL for MocServer mirrors, available in HTTP
  static MIRRORS_HTTPS = [
    "https://alaskybis.unistra.fr/MocServer/query",
    "https://alasky.unistra.fr/MocServer/query"
  ];
  // list of base URL for MocServer mirrors, available in HTTPS
  static _allHiPSes = void 0;
  static _allCatalogHiPSes = void 0;
  static getAllHiPSes() {
    return this._allHiPSes === void 0 && (async () => {
      const A = {
        expr: "dataproduct_type=image||dataproduct_type=cube",
        get: "record",
        fmt: "json",
        fields: "ID,hips_initial_fov,hips_initial_ra,hips_initial_dec,hips_pixel_bitpix,hips_creator,hips_copyright,hips_frame,hips_order,hips_order_min,hips_tile_width,hips_tile_format,hips_pixel_cut,obs_title,obs_description,obs_copyright,obs_regime,hips_data_range,hips_service_url"
      };
      this._allHiPSes = await Z.loadFromMirrors(Lg.MIRRORS_HTTPS, {
        data: A
      }).then((g) => g.json());
    })(), this._allHiPSes;
  }
  static getAllCatalogHiPSes() {
    return this._allCatalogHiPSes === void 0 && (async () => {
      const A = {
        expr: "dataproduct_type=catalog",
        get: "record",
        fmt: "json",
        fields: "ID,hips_copyright,hips_order,hips_order_min,obs_title,obs_description,obs_copyright,obs_regime,cs_service_url,hips_service_url"
      };
      this._allCatalogHiPSes = await Z.loadFromMirrors(Lg.MIRRORS_HTTPS, { data: A }).then((g) => g.json());
    })(), this._allCatalogHiPSes;
  }
}
let Xg = {};
Xg.fetchFromID = async function(B) {
  const A = {
    get: "record",
    fmt: "json",
    ID: "*" + B + "*"
  };
  let g = await Z.loadFromMirrors(Lg.MIRRORS_HTTPS, {
    data: A
  }).then((C) => C.json());
  if (!g || g.length == 0)
    throw "No surveys matching have been found for the id: " + B;
  {
    let C;
    if (g.length > 1) {
      let I = g.find((E) => E.ID === B);
      I ? C = I : (C = g[0], console.warn("Multiple surveys are matching, please choose one. The chosen one is: " + C));
    } else
      C = g[0];
    return C;
  }
};
Xg.fetchFromUrl = async function(B) {
  let A = !1;
  try {
    B = new URL(B);
  } catch {
    try {
      B = Z.getAbsoluteURL(B), B = new URL(B), A = !0;
    } catch (t) {
      throw t;
    }
  }
  const g = B.toString();
  let C = g;
  C.slice(-1) === "/" && (C = C.substr(0, C.length - 1)), C = C + "/properties", A && (C = C + ".txt"), C = Z.getAbsoluteURL(C), C = Z.fixURLForHTTPS(C);
  let I = {};
  return Z.requestCORSIfNotSameOrigin(C) && (I = { mode: "cors" }), await fetch(C, I).then((o) => o.status == 404 ? Promise.reject("Url points to nothing") : o.text()).then((o) => {
    let t = kQ.parseHiPSProperties(o);
    if (t)
      t.hips_service_url = g;
    else
      throw "No surveys matching at this url: " + rootURL;
    return t;
  });
};
Xg.getFasterMirrorUrl = function(B) {
  const A = (I) => {
    I = Z.fixURLForHTTPS(I);
    const E = new AbortController();
    let o = Date.now();
    const t = 2e3, s = setTimeout(() => E.abort(), t);
    return fetch(I + "/properties", { cache: "no-store", signal: E.signal, mode: "cors" }).then((N) => {
      const k = Date.now() - o;
      return clearTimeout(s), { duration: k, baseUrl: I, validRequest: !0 };
    }).catch((N) => ({ duration: t, baseUrl: I, validRequest: !1 }));
  };
  let g = [];
  g.push(A(B.hips_service_url));
  let C = 1;
  for (; B.hasOwnProperty("hips_service_url_" + C.toString()); ) {
    const I = "hips_service_url_" + C.toString();
    let E = B[I];
    g.push(A(E)), C += 1;
  }
  return Promise.all(g).then((I) => {
    let E = I.filter((t) => t.validRequest === !0);
    const o = function(t, s) {
      return t = Math.ceil(t), s = Math.floor(s), Math.floor(Math.random() * (s - t + 1)) + t;
    };
    return E.sort((t, s) => t.duration - s.duration), E.length >= 2 && (E[1].duration - E[0].duration) / E[0].duration < 0.2 ? E[o(0, 1)].baseUrl : E[0].baseUrl;
  }).then((I) => Z.fixURLForHTTPS(I));
};
let mA = {};
mA.tileSize = function(B, A = {}) {
  let g = B && B.tileSize || A.hips_tile_width && +A.hips_tile_width || 512;
  return g & g - 1 !== 0 && (g = 512), g;
};
mA.frame = function(B, A = {}) {
  let g = B && B.cooFrame || A.hips_body && "ICRSd" || A.hips_frame;
  return g == "ICRS" || g == "ICRSd" || g == "equatorial" || g == "j2000" ? g = "ICRS" : g == "galactic" ? g = "GAL" : g === void 0 ? (g = "ICRS", console.warn('No cooframe given. Coordinate systems supported: "ICRS", "ICRSd", "j2000" or "galactic". ICRS is chosen by default')) : (g = "ICRSd", console.warn("Invalid cooframe given: " + cooFrame + '. Coordinate systems supported: "ICRS", "ICRSd", "j2000" or "galactic". ICRS is chosen by default')), g;
};
mA.maxOrder = function(B, A = {}) {
  return B && B.maxOrder || A.hips_order && +A.hips_order;
};
mA.minOrder = function(B, A = {}) {
  return B && B.minOrder || A.hips_order_min && +A.hips_order_min || 0;
};
mA.formats = function(B, A = {}) {
  let g = A.hips_tile_format || "jpeg";
  return g = g.split(" ").map((C) => C.toLowerCase()), g;
};
mA.initialFov = function(B, A = {}) {
  let g = A.hips_initial_fov && +A.hips_initial_fov;
  return g && g < 2777777e-11 && (g = 360), g;
};
mA.skyFraction = function(B, A = {}) {
  return A.moc_sky_fraction && +A.moc_sky_fraction || 0;
};
mA.cutouts = function(B, A = {}) {
  let g = A.hips_pixel_cut && A.hips_pixel_cut.split(" ");
  const C = g && parseFloat(g[0]), I = g && parseFloat(g[1]);
  return [C, I];
};
mA.bitpix = function(B, A = {}) {
  return A.hips_pixel_bitpix && +A.hips_pixel_bitpix;
};
mA.dataproductSubtype = function(B, A = {}) {
  let g = A.dataproduct_subtype || "color";
  return g = g.split(" ").map((C) => C.toLowerCase()), g;
};
mA.isPlanetaryBody = function(B, A = {}) {
  return A.hips_body !== void 0;
};
let EQ = function() {
  function B(g, C, I, E, o) {
    this.view = E, this.wasm = E.wasm, this.added = !1, this.id = g, this.name = C, this.subtype = "survey", this.properties = {}, this.colorCfg = new OI(o);
    let t = this;
    t.query = (async () => {
      let s, M, N, k, U, Y, K, l, d, q, W, T, v, _, a;
      try {
        let z;
        try {
          z = await Xg.fetchFromUrl(I).catch(async (oA) => {
            try {
              return await Xg.fetchFromID(g);
            } catch (hA) {
              throw hA;
            }
          });
        } catch (oA) {
          throw oA;
        }
        if (t.name = t.name || z.obs_title, !z.hips_service_url)
          throw "no valid service URL for retrieving the tiles";
        if (I = Z.fixURLForHTTPS(z.hips_service_url), Xg.getFasterMirrorUrl(z).then((oA) => {
          t.setUrl(oA);
        }), s = mA.maxOrder(o, z), N = mA.tileSize(o, z), k = mA.formats(o, z), d = mA.minOrder(o, z), M = mA.frame(o, z), l = mA.skyFraction(o, z), q = mA.initialFov(o, z), W = +z.hips_initial_ra, T = +z.hips_initial_dec, [U, Y] = mA.cutouts(o, z), K = mA.bitpix(o, z), a = mA.dataproductSubtype(o, z), _ = mA.isPlanetaryBody(o, z), z.hips_body && (v = z.hips_body), z.hips_body !== void 0) {
          if (t.view.options.showFrame) {
            t.view.aladin.setFrame("J2000d");
            let oA = document.querySelectorAll(".aladin-location > .aladin-frameChoice")[0];
            oA.innerHTML = '<option value="' + RA.J2000d.label + '" selected="selected">J2000d</option>';
          }
        } else if (t.view.options.showFrame) {
          const oA = RA.fromString(t.view.options.cooFrame, RA.J2000);
          let hA = document.querySelectorAll(".aladin-location > .aladin-frameChoice")[0];
          hA.innerHTML = '<option value="' + RA.J2000.label + '" ' + (oA == RA.J2000 ? 'selected="selected"' : "") + '>J2000</option><option value="' + RA.J2000d.label + '" ' + (oA == RA.J2000d ? 'selected="selected"' : "") + '>J2000d</option><option value="' + RA.GAL.label + '" ' + (oA == RA.GAL ? 'selected="selected"' : "") + ">GAL</option>";
        }
      } catch (z) {
        throw z;
      }
      t.properties = {
        url: I,
        maxOrder: s,
        frame: M,
        tileSize: N,
        formats: k,
        minCutout: U,
        maxCutout: Y,
        bitpix: K,
        skyFraction: l,
        minOrder: d,
        hipsInitialFov: q,
        hipsInitialRa: W,
        hipsInitialDec: T,
        dataproductSubtype: a,
        isPlanetaryBody: _,
        hipsBody: v
      };
      let iA = !1;
      o && o.longitudeReversed === !0 && (iA = !0), t.properties.hipsBody && (iA = !0), t.longitudeReversed = iA;
      let IA = o && o.imgFormat;
      if (IA) {
        if (IA = IA.toLowerCase(), IA === "jpg" && (IA = "jpeg"), IA === "fits" && k.indexOf("fits") < 0)
          throw t.name + " does not provide fits tiles";
        if (IA === "webp" && k.indexOf("webp") < 0)
          throw t.name + " does not provide webp tiles";
        if (IA === "png" && k.indexOf("png") < 0)
          throw t.name + " does not provide png tiles";
        if (IA === "jpeg" && k.indexOf("jpeg") < 0)
          throw t.name + " does not provide jpeg tiles";
      } else if (k.indexOf("png") >= 0)
        IA = "png";
      else if (k.indexOf("webp") >= 0)
        IA = "webp";
      else if (k.indexOf("jpeg") >= 0)
        IA = "jpeg";
      else if (k.indexOf("fits") >= 0)
        IA = "fits";
      else
        throw "Unsupported format(s) found in the properties: " + k;
      if (t.imgFormat = IA, IA === "fits") {
        const z = o && o.minCut || t.properties.minCutout || 0, oA = o && o.maxCut || t.properties.maxCutout || 1;
        this.colorCfg.setCuts(z, oA);
      }
      return Fg.update(t), t;
    })();
  }
  B.prototype.isReady = function() {
    return this.added;
  }, B.prototype.setUrl = function(g) {
    this.properties.url !== g && (console.info("Change url of ", this.id, " from ", this.properties.url, " to ", g), this.added && this.wasm.setHiPSUrl(this.properties.url, g), this.properties.url = g);
  }, B.prototype.isPlanetaryBody = function() {
    return this.properties.isPlanetaryBody;
  }, B.prototype.setImageFormat = function(g) {
    let C = this;
    C.query.then(() => {
      A(C, () => {
        let I = g.toLowerCase();
        if (I !== "fits" && I !== "png" && I !== "jpg" && I !== "jpeg" && I !== "webp")
          throw 'Formats must lie in ["fits", "png", "jpg", "webp"]';
        if (I === "jpg" && (I = "jpeg"), C.imgFormat === I)
          return;
        const E = C.properties.formats;
        if (I === "fits" && E.indexOf("fits") < 0)
          throw C.id + " does not provide fits tiles";
        if (I === "webp" && E.indexOf("webp") < 0)
          throw C.id + " does not provide webp tiles";
        if (I === "png" && E.indexOf("png") < 0)
          throw C.id + " does not provide png tiles";
        if (I === "jpeg" && E.indexOf("jpeg") < 0)
          throw C.id + " does not provide jpeg tiles";
        C.imgFormat = I;
      });
    });
  }, B.prototype.setOpacity = function(g) {
    let C = this;
    A(C, () => {
      C.colorCfg.setOpacity(g);
    });
  }, B.prototype.setBlendingConfig = function(g = !1) {
    A(this, () => {
      this.colorCfg.setBlendingConfig(g);
    });
  }, B.prototype.setColormap = function(g, C) {
    A(this, () => {
      this.colorCfg.setColormap(g, C);
    });
  }, B.prototype.setCuts = function(g, C) {
    A(this, () => {
      this.colorCfg.setCuts(g, C);
    });
  }, B.prototype.setGamma = function(g) {
    A(this, () => {
      this.colorCfg.setGamma(g);
    });
  }, B.prototype.setSaturation = function(g) {
    A(this, () => {
      this.colorCfg.setSaturation(g);
    });
  }, B.prototype.setBrightness = function(g) {
    A(this, () => {
      this.colorCfg.setBrightness(g);
    });
  }, B.prototype.setContrast = function(g) {
    A(this, () => {
      this.colorCfg.setContrast(g);
    });
  }, B.prototype.metadata = function() {
    return {
      ...this.colorCfg.get(),
      longitudeReversed: this.longitudeReversed,
      imgFormat: this.imgFormat
    };
  };
  var A = function(g, C = void 0) {
    C && C();
    try {
      if (g.added) {
        const I = g.metadata();
        g.wasm.setImageMetadata(g.layer, I), CA.HIPS_LAYER_CHANGED.dispatchedTo(g.view.aladinDiv, { layer: g });
      }
    } catch (I) {
      console.error(I);
    }
  };
  return B.prototype.add = function(g) {
    return this.layer = g, this.wasm.addImageSurvey({
      layer: this.layer,
      properties: this.properties,
      meta: this.metadata()
    }), this.added = !0, Promise.resolve(this);
  }, B.prototype.toggle = function() {
    this.colorCfg.getOpacity() != 0 ? this.colorCfg.setOpacity(0) : this.colorCfg.setOpacity(this.prevOpacity);
  }, B.prototype.setAlpha = B.prototype.setOpacity, B.prototype.setColorCfg = function(g) {
    A(this, () => {
      this.colorCfg = g;
    });
  }, B.prototype.getColorCfg = function() {
    return this.colorCfg;
  }, B.prototype.getOpacity = function() {
    return this.colorCfg.getOpacity();
  }, B.prototype.getAlpha = B.prototype.getOpacity, B.prototype.readPixel = function(g, C) {
    return this.wasm.readPixel(g, C, this.layer);
  }, B.DEFAULT_SURVEY_ID = "P/DSS2/color", B;
}(), YE = function() {
  let B = {};
  return B.GALACTIC_TO_J2000 = [
    -0.0548755604024359,
    0.4941094279435681,
    -0.867666148981161,
    -0.8734370902479237,
    -0.4448296299195045,
    -0.1980763734646737,
    -0.4838350155267381,
    0.7469822444763707,
    0.4559837762325372
  ], B.J2000_TO_GALACTIC = [
    -0.0548755604024359,
    -0.873437090247923,
    -0.4838350155267381,
    0.4941094279435681,
    -0.4448296299195045,
    0.7469822444763707,
    -0.867666148981161,
    -0.1980763734646737,
    0.4559837762325372
  ], B.Transform = function(A, g) {
    A[0] = A[0] * Math.PI / 180, A[1] = A[1] * Math.PI / 180;
    var C = new Array(
      Math.cos(A[0]) * Math.cos(A[1]),
      Math.sin(A[0]) * Math.cos(A[1]),
      Math.sin(A[1])
    ), I = new Array(
      C[0] * g[0] + C[1] * g[1] + C[2] * g[2],
      C[0] * g[3] + C[1] * g[4] + C[2] * g[5],
      C[0] * g[6] + C[1] * g[7] + C[2] * g[8]
    ), E = Math.sqrt(I[0] * I[0] + I[1] * I[1] + I[2] * I[2]), o = new Array(0, 0);
    o[1] = Math.asin(I[2] / E);
    var t = I[0] / E / Math.cos(o[1]), s = I[1] / E / Math.cos(o[1]);
    return o[0] = Math.atan2(s, t), o[0] < 0 && (o[0] = o[0] + 2 * Math.PI), o[0] = o[0] * 180 / Math.PI, o[1] = o[1] * 180 / Math.PI, o;
  }, B.GalacticToJ2000 = function(A) {
    return B.Transform(A, B.GALACTIC_TO_J2000);
  }, B.J2000ToGalactic = function(A) {
    return B.Transform(A, B.J2000_TO_GALACTIC);
  }, B;
}();
const lE = function() {
  return function(A) {
    const g = document.createElement("div");
    g.classList.add("aladin-logo-container");
    const C = document.createElement("a");
    C.href = "https://aladin.cds.unistra.fr/", C.title = "Powered by Aladin Lite", C.target = "_blank";
    const I = document.createElement("div");
    I.classList.add("aladin-logo"), C.appendChild(I), g.appendChild(C), A.appendChild(g);
  };
}();
class JE {
  // Constructor
  constructor(A, g) {
    this.aladin = g, this.mainDiv = document.createElement("div"), this.mainDiv.classList.add("aladin-projection-select"), A.appendChild(this.mainDiv), this._createComponent(), this._addListeners();
  }
  _createComponent() {
    O(this.mainDiv).append('<select title="Projection" class="aladin-selector"></select>'), this.selectProjection = O(this.mainDiv).find("select"), this.selectProjection.empty(), pC.forEach((g) => {
      this.selectProjection.append(O("<option />").val(g).text(g));
    });
    let A = this;
    this.selectProjection.change(function() {
      A.aladin.setProjection(O(this).val());
    });
  }
  _addListeners() {
    const A = this;
    CA.PROJECTION_CHANGED.listenedBy(this.aladin.aladinDiv, function(g) {
      A.selectProjection.val(g.detail.projection);
    });
  }
  show() {
    this.mainDiv.style.display = "block";
  }
  hide() {
    this.mainDiv.style.display = "none";
  }
}
var nQ = { exports: {} };
(function(B, A) {
  (function(g, C) {
    B.exports = C();
  })(wQ, function() {
    function g(C) {
      var I = document, E = C.container || I.createElement("div"), o = E.style, t = navigator.userAgent, s = ~t.indexOf("Firefox") && ~t.indexOf("Mobile"), M = C.debounceWaitMs || 0, N = C.preventSubmit || !1, k = C.disableAutoSelect || !1, U = s ? "input" : "keyup", Y = [], K = "", l = 2, d = C.showOnFocus, q, W = 0, T;
      if (C.minLength !== void 0 && (l = C.minLength), !C.input)
        throw new Error("input undefined");
      var v = C.input;
      E.className = "autocomplete " + (C.className || ""), o.position = "absolute";
      function _() {
        var sA = E.parentNode;
        sA && sA.removeChild(E);
      }
      function a() {
        T && window.clearTimeout(T);
      }
      function iA() {
        E.parentNode || I.body.appendChild(E);
      }
      function IA() {
        return !!E.parentNode;
      }
      function z() {
        W++, Y = [], K = "", q = void 0, _();
      }
      function oA() {
        if (!IA())
          return;
        o.height = "auto", o.width = v.offsetWidth + "px";
        var sA = 0, KA;
        function xA() {
          var vA = I.documentElement, ZA = vA.clientTop || I.body.clientTop || 0, HA = vA.clientLeft || I.body.clientLeft || 0, Rg = window.pageYOffset || vA.scrollTop, Ig = window.pageXOffset || vA.scrollLeft;
          KA = v.getBoundingClientRect();
          var qg = KA.top + v.offsetHeight + Rg - ZA, ag = KA.left + Ig - HA;
          o.top = qg + "px", o.left = ag + "px", sA = window.innerHeight - (KA.top + v.offsetHeight), sA < 0 && (sA = 0), o.top = qg + "px", o.bottom = "", o.left = ag + "px", o.maxHeight = sA + "px";
        }
        xA(), xA(), C.customize && KA && C.customize(v, KA, E, sA);
      }
      function hA() {
        for (; E.firstChild; )
          E.removeChild(E.firstChild);
        var sA = function(HA, Rg) {
          var Ig = I.createElement("div");
          return Ig.textContent = HA.label || "", Ig;
        };
        C.render && (sA = C.render);
        var KA = function(HA, Rg) {
          var Ig = I.createElement("div");
          return Ig.textContent = HA, Ig;
        };
        C.renderGroup && (KA = C.renderGroup);
        var xA = I.createDocumentFragment(), vA = "#9?$";
        if (Y.forEach(function(HA) {
          if (HA.group && HA.group !== vA) {
            vA = HA.group;
            var Rg = KA(HA.group, K);
            Rg && (Rg.className += " group", xA.appendChild(Rg));
          }
          var Ig = sA(HA, K);
          Ig && (Ig.addEventListener("click", function(qg) {
            C.onSelect(HA, v), z(), qg.preventDefault(), qg.stopPropagation();
          }), HA === q && (Ig.className += " selected"), xA.appendChild(Ig));
        }), E.appendChild(xA), Y.length < 1)
          if (C.emptyMsg) {
            var ZA = I.createElement("div");
            ZA.className = "empty", ZA.textContent = C.emptyMsg, E.appendChild(ZA);
          } else {
            z();
            return;
          }
        iA(), oA(), og();
      }
      function yA() {
        IA() && hA();
      }
      function rA() {
        yA();
      }
      function gA(sA) {
        sA.target !== E ? yA() : sA.preventDefault();
      }
      function nA(sA) {
        for (var KA = sA.which || sA.keyCode || 0, xA = C.keysToIgnore || [
          38,
          13,
          27,
          39,
          37,
          16,
          17,
          18,
          20,
          91,
          9
          /* Tab */
        ], vA = 0, ZA = xA; vA < ZA.length; vA++) {
          var HA = ZA[vA];
          if (KA === HA)
            return;
        }
        KA >= 112 && KA <= 123 && !C.keysToIgnore || KA === 40 && IA() || Dg(
          0
          /* Keyboard */
        );
      }
      function og() {
        var sA = E.getElementsByClassName("selected");
        if (sA.length > 0) {
          var KA = sA[0], xA = KA.previousElementSibling;
          if (xA && xA.className.indexOf("group") !== -1 && !xA.previousElementSibling && (KA = xA), KA.offsetTop < E.scrollTop)
            E.scrollTop = KA.offsetTop;
          else {
            var vA = KA.offsetTop + KA.offsetHeight, ZA = E.scrollTop + E.offsetHeight;
            vA > ZA && (E.scrollTop += vA - ZA);
          }
        }
      }
      function PA() {
        if (Y.length < 1)
          q = void 0;
        else if (q === Y[0])
          q = Y[Y.length - 1];
        else
          for (var sA = Y.length - 1; sA > 0; sA--)
            if (q === Y[sA] || sA === 1) {
              q = Y[sA - 1];
              break;
            }
      }
      function gg() {
        if (Y.length < 1 && (q = void 0), !q || q === Y[Y.length - 1]) {
          q = Y[0];
          return;
        }
        for (var sA = 0; sA < Y.length - 1; sA++)
          if (q === Y[sA]) {
            q = Y[sA + 1];
            break;
          }
      }
      function XA(sA) {
        var KA = sA.which || sA.keyCode || 0;
        if (KA === 38 || KA === 40 || KA === 27) {
          var xA = IA();
          if (KA === 27)
            z();
          else {
            if (!xA || Y.length < 1)
              return;
            KA === 38 ? PA() : gg(), hA();
          }
          sA.preventDefault(), xA && sA.stopPropagation();
          return;
        }
        KA === 13 && (q && (C.onSelect(q, v), z()), N && sA.preventDefault());
      }
      function OA() {
        d && Dg(
          1
          /* Focus */
        );
      }
      function Dg(sA) {
        var KA = ++W, xA = v.value, vA = v.selectionStart || 0;
        xA.length >= l || sA === 1 ? (a(), T = window.setTimeout(function() {
          C.fetch(xA, function(ZA) {
            W === KA && ZA && (Y = ZA, K = xA, q = Y.length < 1 || k ? void 0 : Y[0], hA());
          }, sA, vA);
        }, sA === 0 ? M : 0)) : z();
      }
      function RI() {
        setTimeout(function() {
          I.activeElement !== v && z();
        }, 200);
      }
      E.addEventListener("mousedown", function(sA) {
        sA.stopPropagation(), sA.preventDefault();
      }), E.addEventListener("focus", function() {
        return v.focus();
      });
      function Jg() {
        v.removeEventListener("focus", OA), v.removeEventListener("keydown", XA), v.removeEventListener(U, nA), v.removeEventListener("blur", RI), window.removeEventListener("resize", rA), I.removeEventListener("scroll", gA, !0), a(), z();
      }
      return v.addEventListener("keydown", XA), v.addEventListener(U, nA), v.addEventListener("blur", RI), v.addEventListener("focus", OA), window.addEventListener("resize", rA), I.addEventListener("scroll", gA, !0), {
        destroy: Jg
      };
    }
    return g;
  });
})(nQ);
var UE = nQ.exports;
const FQ = /* @__PURE__ */ tQ(UE);
class SE {
  constructor(A, g, C) {
    this.parentDiv = A, this.aladin = g, this.fnIdSelected = C, this._createComponent();
  }
  _createComponent() {
    const A = this;
    this.mainDiv = document.createElement("div"), this.mainDiv.classList.add("aladin-dialog", "aladin-cb-list"), this.mainDiv.style.display = "block";
    const g = "autocomplete-" + Z.uuidv4();
    this.mainDiv.insertAdjacentHTML(
      "afterbegin",
      '<a class="aladin-closeBtn">&times;</a><div class="aladin-box-title">Select Catalogue:</div><div class="aladin-box-content"><div class="aladin-label" for="' + g + '">By ID, title, keyword</div><input class="aladin-input" style="width:100%;" name="' + g + '" id="' + g + '" type="text" placeholder="Type keyword or VOTable URL" /><div class="aladin-row cone-search" style="font-size: 12px;"><div><input class="aladin-input" type="number" value="1.0" style="width: 4em;" maxlength="5" size="5"> <select class="aladin-selector" style="padding: 4px 0!important;"><option>deg<option>arcmin<option>arcsec</select> around view center</div><div>Limit to <input class="aladin-input" type="number" min="1" max="10000" value="1000" style="width: 5em;"> sources</div></div><div class="aladin-row"><div class="aladin-col"><div><button class="aladin-btn">Load cone</button></div></div><div class="hips aladin-col"><button class="aladin-btn">Load catalogue HiPS</button></div></div><div class="aladin-box-separator"></div><div class="aladin-label" for="' + g + '">By VOTable URL</div><input class="aladin-input" style="width:100%;" name="' + g + '" id="' + g + '" type="text" placeholder="Enter VOTable URL" /><div class="votable"><button class="aladin-btn">Load VOTable</button></div></div>'
    ), this.parentDiv.appendChild(this.mainDiv), this.idInput = A.mainDiv.querySelectorAll("input")[0], this.votInput = A.mainDiv.querySelectorAll("input")[3], O(this.idInput).on("change", function() {
      A.idInput.blur();
    }), O(this.votInput).on("change", function() {
      A.votInput.blur();
    });
    let [C, I, E] = this.mainDiv.querySelectorAll(".aladin-btn");
    this.divCS = this.mainDiv.querySelector(".cone-search"), this.divLoadHiPS = this.mainDiv.querySelector(".hips"), this.divLoadHiPS.style.display = "none", O(C).click(function() {
      const s = parseFloat(A.divCS.querySelector("div:nth-child(1) > input").value), M = A.divCS.querySelector("div:nth-child(1) > select").value;
      let N = s;
      M == "arcmin" ? N /= 60 : M == "arcsec" && (N /= 3600);
      const k = parseInt(A.divCS.querySelector("div:nth-child(2) > input").value), U = A.selectedItem.cs_service_url, [Y, K] = A.aladin.getRaDec();
      A.fnIdSelected && A.fnIdSelected("coneSearch", { id: A.idInput.value, baseURL: U, limit: k, radiusDeg: N, ra: Y, dec: K }), A.idInput.value = null;
    }), O(I).click(function() {
      A.fnIdSelected && A.fnIdSelected("hips", { id: A.idInput.value, hipsURL: A.selectedItem.hips_service_url }), A.idInput.value = null;
    }), O(E).click(function() {
      A.fnIdSelected && A.fnIdSelected("votable", { url: A.votInput.value }), A.votInput.value = null;
    });
    let o = document.getElementById(g);
    Lg.getAllCatalogHiPSes(), FQ({
      input: o,
      minLength: 3,
      fetch: function(s, M) {
        s = s.toLowerCase();
        const N = function(Y) {
          const K = Y.ID, l = Y.obs_title || "", d = Y.obs_description || "";
          return K.toLowerCase().includes(s) || l.toLowerCase().includes(s) || d.toLowerCase().includes(s);
        }, k = Lg.getAllCatalogHiPSes().filter(N);
        k.sort(function(Y, K) {
          let l = 0, d = 0;
          return Y.ID.toLowerCase().includes(s) && (l += 100), K.ID.toLowerCase().includes(s) && (d += 100), Y.obs_title.toLowerCase().includes(s) && (l += 50), K.obs_title.toLowerCase().includes(s) && (d += 50), Y.obs_description.toLowerCase().includes(s) && (l += 10), K.obs_description.toLowerCase().includes(s) && (d += 10), Y.hips_service_url && (l += 20), K.hips_service_url && (d += 20), l > d ? -1 : d > l ? 1 : 0;
        });
        const U = k.slice(0, 50);
        M(U);
      },
      onSelect: function(s) {
        s.cs_service_url ? O(A.divCS).show() : O(A.divCS).hide(), s.hips_service_url ? O(A.divLoadHiPS).show() : O(A.divLoadHiPS).hide(), o.value = s.ID, A.selectedItem = s, o.blur();
      },
      // attach container to AL div if needed (to prevent it from being hidden in full screen mode)
      customize: function(s, M, N, k) {
        A.aladin.fullScreenBtn.hasClass("aladin-restore") && A.parentDiv.appendChild(N);
      },
      render: function(s, M) {
        const N = document.createElement("div");
        return N.innerHTML = (s.obs_title || "") + ' - <span style="color: #ae8de1">' + s.ID + "</span>", N;
      }
    });
    let [t] = this.mainDiv.querySelectorAll(".aladin-closeBtn");
    O(t).click(function() {
      A.hide();
    });
  }
  show() {
    this.mainDiv.style.display = "block";
  }
  hide() {
    this.mainDiv.style.display = "none";
  }
}
class LE {
  constructor(A, g, C) {
    this.parentDiv = A, this.fnIdSelected = g, this.aladin = C, this._createComponent();
  }
  _createComponent() {
    const A = this;
    this.mainDiv = document.createElement("div"), this.mainDiv.style.display = "block", this.mainDiv.classList.add("aladin-dialog", "aladin-layerBox", "aladin-cb-list");
    const g = "autocomplete-" + Z.uuidv4();
    this.mainDiv.insertAdjacentHTML(
      "afterbegin",
      '<a class="aladin-closeBtn">&times;</a><div class="aladin-box-title">Select image HiPS:</div><div class="aladin-box-content"><div class="aladin-label" for="' + g + '">By ID, title, keyword or URL</div><input class="aladin-input" style="width:100%" name="' + g + '" id="' + g + '" type="text" placeholder="Type ID, title, keyword or URL" /><br><div><button class="aladin-btn">Select HiPS</button><button class="aladin-btn">Load coverage</button></div></div>'
    ), this.parentDiv.appendChild(this.mainDiv);
    const C = document.getElementById(g);
    O(C).on("change", function() {
      C.blur();
    }), Lg.getAllHiPSes(), FQ({
      input: C,
      fetch: function(t, s) {
        t = t.toLowerCase();
        const M = Lg.getAllHiPSes().filter((k) => k.ID.toLowerCase().includes(t) || k.obs_title.toLowerCase().includes(t));
        M.sort(function(k, U) {
          let Y = 0, K = 0;
          return k.ID.toLowerCase().includes(t) && (Y += 100), U.ID.toLowerCase().includes(t) && (K += 100), k.obs_title.toLowerCase().includes(t) && (Y += 50), U.obs_title.toLowerCase().includes(t) && (K += 50), k.obs_description && k.obs_description.toLowerCase().includes(t) && (Y += 10), U.obs_description && U.obs_description.toLowerCase().includes(t) && (K += 10), Y > K ? -1 : K > Y ? 1 : 0;
        });
        const N = M.slice(0, 50);
        s(N);
      },
      onSelect: function(t) {
        A.selectedItem = t, C.value = t.ID, A.fnIdSelected && A.fnIdSelected(t.ID), C.blur();
      },
      // attach container to AL div if needed (to prevent it from being hidden in full screen mode)
      customize: function(t, s, M, N) {
        A.aladin.fullScreenBtn.hasClass("aladin-restore") && A.parentDiv.appendChild(M);
      },
      render: function(t, s) {
        const M = document.createElement("div");
        return M.innerHTML = t.obs_title + ' - <span style="color: #ae8de1">' + t.ID + "</span>", M;
      }
    });
    let [I, E] = this.mainDiv.querySelectorAll(".aladin-btn"), [o] = this.mainDiv.querySelectorAll(".aladin-closeBtn");
    O(o).click(function() {
      A.hide();
    }), O(I).click(function() {
      let t = A.mainDiv.querySelectorAll("input")[0];
      t && A.fnIdSelected && A.fnIdSelected(t.value), t.value = "";
    }), O(E).on("click", function() {
      let t, s = A.mainDiv.querySelectorAll("input")[0];
      s.value.startsWith("http") ? t = s.value + "/Moc.fits" : t = A.selectedItem.hips_service_url + "/Moc.fits", t = Z.fixURLForHTTPS(t);
      const M = zA.MOCFromURL(t, { lineWidth: 5, opacity: 0.3 });
      A.aladin.addMOC(M);
    });
  }
  show() {
    this.mainDiv.style.display = "block";
  }
  hide() {
    this.mainDiv.style.display = "none";
  }
}
class iQ {
  // Constructor
  constructor(A, g) {
    if (this.aladin = A, this.layer = g, this.hidden = !1, this.lastOpacity = 1, this.headerDiv = O(
      '<div class="aladin-layer"><div class="aladin-layer-header" style="border-radius: 4px"><button class="aladin-btn-small aladin-indicatorBtn right-triangle" title="Open the color panel"></button><select class="aladin-selector aladin-layerSelection"></select><button class="aladin-btn-small aladin-layer-hide" type="button" title="Hide this layer">👁️</button><button class="aladin-btn-small aladin-HiPSSelector" type="button" title="Search for a specific HiPS">🔍</button><button class="aladin-btn-small aladin-delete-layer" type="button" title="Delete this layer">❌</button></div></div>'
    ), this.layer.subtype === "fits" && this.headerDiv[0].querySelector(".aladin-layerSelection").after(O('<button class="aladin-btn-small aladin-layer-focuson" type="button" title="Focus on this layer">🎯</button>')[0]), this.layer.layer === "base") {
      let E = this.headerDiv[0].querySelector(".aladin-delete-layer");
      E.disabled = !0, E.style.backgroundColor = "lightgray", E.style.borderColor = "gray", E.style.color = "transparent", E.style.textShadow = "0 0 0 gray";
    }
    let C = "";
    for (const E of this.aladin.wasm.getAvailableColormapList())
      C += "<option>" + E + "</option>";
    this.cmap = "native", this.color = "#ff0000", this.mainDiv = O('<div class="aladin-frame" style="display:none; padding: 0px 4px"><div class="aladin-options">  <div class="row"><div class="col-label">Colormap</div><div class="col-input"><select class="aladin-selector colormap-selector">' + C + '</select></div></div>  <label><div class="row"><div class="col-label">Reverse</div><div class="col-input"><input type="checkbox" class="reversed aladin-input"></div></div></label>  <div class="row"><div class="col-label"><label>Stretch</label></div><div class="col-input"><select class="aladin-selector stretch"><option>pow2</option><option selected>linear</option><option>sqrt</option><option>asinh</option><option>log</option></select></div></div>  <div class="row"><div class="col-label"><label>Format</label></div><div class="col-input"><select class="aladin-selector format"></select></div></div>  <div class="row"><div class="col-label"><label>Min cut</label></div><div class="col-input"><input type="number" class="aladin-input min-cut"></div></div>  <div class="row"><div class="col-label"><label>Max cut</label></div><div class="col-input"><input type="number" class="aladin-input max-cut"></div></div>  <div class="row"><div class="col-label"><label>Gamma</label></div><div class="col-input"><input class="aladin-input gamma" type="number" value="1.0" min="0.1" max="10.0" step="0.01"></div></div>  <div class="row"><div class="col-label"><label>Color Sat.</label></div><div class="col-input"><input class="aladin-input saturation" type="range" value="0.0" min="-1.0" max="1.0" step="0.01"></div></div>  <div class="row"><div class="col-label"><label>Contrast</label></div><div class="col-input"><input class="aladin-input contrast" type="range" value="0.0" min="-1.0" max="1.0" step="0.01"></div></div>  <div class="row"><div class="col-label"><label>Brightness</label></div><div class="col-input"><input class="aladin-input brightness" type="range" value="0.0" min="-1.0" max="1.0" step="0.01"></div></div>  <div class="row"><div class="col-label"><label>Blending mode</label></div><div class="col-input"><select class="aladin-selector blending"><option>additive</option><option selected>default</option></select></div></div>  <div class="row"><div class="col-label"><label>Opacity</label></div><div class="col-input"><input class="aladin-input opacity" type="range" min="0" max="1" step="0.01"></div></div></div> </div>'), this._addListeners(), this._updateHiPSLayerOptions();
    let I = this;
    this.layerChangedListener = function(E) {
      const o = E.detail.layer;
      o.layer === I.layer.layer && (I.layer = o, I._updateHiPSLayerOptions()), I._updateLayersDropdownList();
    }, CA.HIPS_LAYER_CHANGED.listenedBy(this.aladin.aladinDiv, this.layerChangedListener);
  }
  destroy() {
    CA.HIPS_LAYER_CHANGED.remove(this.aladin.aladinDiv, this.layerChangedListener);
  }
  _addListeners() {
    const A = this, g = this.headerDiv.find(".aladin-indicatorBtn");
    g.off("click"), g.on("click", function() {
      g.hasClass("right-triangle") ? (g.removeClass("right-triangle"), g.addClass("down-triangle"), A.mainDiv.slideDown(300)) : (g.removeClass("down-triangle"), g.addClass("right-triangle"), A.mainDiv.slideUp(300));
    }), this.headerDiv.off("click"), this.headerDiv.on("click", () => {
      A.aladin.aladinDiv.dispatchEvent(new CustomEvent("select-layer", {
        detail: A
      }));
    }), A._updateLayersDropdownList();
    const C = this.headerDiv.find(".aladin-layerSelection");
    C.off("change"), C.on("change", (v) => {
      let _ = Fg.LAYERS[C[0].selectedIndex], a;
      _.subtype === "fits" ? a = A.aladin.createImageFITS(
        _.url,
        _.name,
        _.options
      ) : a = A.aladin.createImageSurvey(
        _.id,
        _.name,
        _.url,
        void 0,
        _.maxOrder,
        _.options
      ), A.hidden && a.setAlpha(0), A.aladin.setOverlayImageLayer(a, A.layer.layer);
    });
    const I = this.headerDiv.find(".aladin-HiPSSelector");
    I.off("click"), I.on("click", function() {
      A.hipsSelector || (A.hipsSelector = new LE(A.aladin.aladinDiv, (v) => {
        const _ = A.layer.layer;
        A.aladin.setOverlayImageLayer(v, _);
      }, A.aladin)), A.hipsSelector.show();
    });
    const E = this.headerDiv.find(".aladin-delete-layer");
    E.off("click"), E.on("click", function() {
      const v = new CustomEvent("remove-layer", {
        detail: A.layer.layer
      });
      A.aladin.aladinDiv.dispatchEvent(v);
    });
    const o = this.headerDiv.find(".aladin-layer-hide");
    o.off("click"), o.on("click", function() {
      A.hidden = !A.hidden;
      let v = A.mainDiv.find(".opacity").eq(0), _ = 0;
      A.hidden ? (A.lastOpacity = A.layer.getOpacity(), o.text("")) : (_ = A.lastOpacity, o.text("👁️")), v.val(_), v.get(0).disabled = A.hidden, A.layer.setOpacity(_);
    });
    const t = this.headerDiv.find(".aladin-layer-focuson");
    t && (t.off("click"), t.on("click", function() {
      A.layer.focusOn();
    }));
    const s = this.mainDiv.find(".blending").eq(0);
    s.off("change"), s.change(function() {
      let v = s.val();
      A.layer.setBlendingConfig(v === "additive");
    });
    const M = this.mainDiv.find(".format").eq(0), N = this.mainDiv.find(".min-cut").eq(0), k = this.mainDiv.find(".max-cut").eq(0);
    M.off("change"), M.on("change", function() {
      const v = M.val();
      A.layer.setImageFormat(v);
      let _ = 0, a = 1;
      v === "fits" && (_ = A.layer.properties.minCutout, a = A.layer.properties.maxCutout), A.layer.setCuts(_, a), N.val(parseFloat(_.toFixed(5))), k.val(parseFloat(a.toFixed(5)));
    }), N.off("change blur"), k.off("change blur"), N.add(k).on("change blur", function(v) {
      let _ = parseFloat(N.val()), a = parseFloat(k.val());
      isNaN(_) || isNaN(a) || A.layer.setCuts(_, a);
    });
    const U = this.mainDiv.find(".colormap-selector").eq(0), Y = this.mainDiv.find(".stretch").eq(0), K = this.mainDiv.find(".reversed").eq(0);
    K.off("change"), U.off("change"), Y.off("change"), U.add(K).add(Y).change(function() {
      const v = Y.val(), _ = K[0].checked, a = U.val();
      A.layer.setColormap(a, { reversed: _, stretch: v });
    });
    const l = A.mainDiv.find(".opacity").eq(0);
    l.off("input"), l.on("input", function() {
      const v = +l.val();
      A.layer.setOpacity(v);
    });
    const d = A.mainDiv.find(".gamma").eq(0);
    d.off("change blur"), d.on("change blur", function() {
      const v = parseFloat(d.val()) || 1;
      A.layer.setGamma(v);
      const _ = A.layer.getColorCfg().getGamma();
      v !== _ && d.val(_);
    });
    const q = A.mainDiv.find(".saturation").eq(0);
    q.off("input"), q.on("input", function(v) {
      const _ = parseFloat(q.val()) || 0;
      A.layer.setSaturation(_);
      const a = A.layer.getColorCfg().getSaturation();
      _ !== a && q.val(a);
    });
    const W = A.mainDiv.find(".contrast").eq(0);
    W.off("input"), W.on("input", function(v) {
      const _ = parseFloat(W.val()) || 0;
      A.layer.setContrast(_);
      const a = A.layer.getColorCfg().getContrast();
      _ !== a && W.val(a);
    });
    const T = A.mainDiv.find(".brightness").eq(0);
    T.off("input"), T.on("input", function(v) {
      const _ = parseFloat(T.val()) || 0;
      A.layer.setBrightness(_);
      const a = A.layer.getColorCfg().getBrightness();
      _ !== a && T.val(a);
    });
  }
  _updateHiPSLayerOptions() {
    const A = this.mainDiv.find(".row").eq(0), g = this.mainDiv.find(".row").eq(1), C = this.mainDiv.find(".row").eq(2);
    this.mainDiv.find(".row").eq(3);
    const I = this.mainDiv.find(".row").eq(4), E = this.mainDiv.find(".row").eq(5), o = this.mainDiv.find(".reversed").eq(0), t = this.mainDiv.find(".colormap-selector").eq(0), s = this.mainDiv.find(".stretch").eq(0), M = this.mainDiv.find(".format").eq(0), N = this.mainDiv.find(".opacity").eq(0), k = this.mainDiv.find(".gamma").eq(0), U = this.mainDiv.find(".contrast").eq(0), Y = this.mainDiv.find(".brightness").eq(0), K = this.mainDiv.find(".saturation").eq(0), l = this.mainDiv.find(".blending").eq(0), d = this.mainDiv.find(".min-cut").eq(0), q = this.mainDiv.find(".max-cut").eq(0);
    M.empty(), this.layer.properties.formats.forEach((rA) => {
      M.append(O("<option>", {
        value: rA,
        text: rA
      }));
    });
    const W = this.layer.getColorCfg(), T = W.colormap, v = W.reversed, _ = W.stretch, a = this.layer.imgFormat;
    M.val(a);
    const iA = W.getBlendingConfig();
    l.val(iA ? "additive" : "default"), A[0].style.display = "flex", g[0].style.display = "flex", C[0].style.display = "flex", W.minCut ? parseFloat(d.val()) != W.minCut && d.val(parseFloat(W.minCut.toFixed(5))) : d.val(0), I[0].style.display = "flex", W.maxCut ? parseFloat(q.val()) != W.maxCut && q.val(parseFloat(W.maxCut.toFixed(5))) : q.val(0), E[0].style.display = "flex";
    const IA = W.getOpacity();
    N.val(IA);
    const z = W.getGamma();
    k.val(z);
    const oA = W.getSaturation();
    K.val(oA);
    const hA = W.getBrightness();
    Y.val(hA);
    const yA = W.getContrast();
    U.val(yA), t.val(T), this.cmap = T, o.prop("checked", v), s.val(_);
  }
  _updateLayersDropdownList() {
    let A = this.headerDiv.find(".aladin-layerSelection"), g = Fg.LAYERS.sort(function(C, I) {
      return C.order ? C.maxOrder && C.maxOrder > I.maxOrder ? 1 : -1 : C.name > I.name ? 1 : -1;
    });
    A.empty(), this.layer && g.forEach((C) => {
      const I = this.layer.id.endsWith(C.id);
      A.append(O("<option />").attr("selected", I).val(C.id).text(C.name));
    });
  }
  attachTo(A) {
    this.headerDiv.append(this.mainDiv), A.append(this.headerDiv), this._addListeners();
  }
  show() {
    this.mainDiv.style.display = "block";
  }
  hide() {
    this.headerDiv.style.display = "none", this.mainDiv.style.display = "none";
  }
}
class qE {
  // Constructor
  constructor(A, g, C) {
    this.aladin = g, this.view = C, this.mainDiv = document.createElement("div"), this.mainDiv.style.display = "none", this.mainDiv.classList.add("aladin-box", "aladin-layerBox", "aladin-cb-list"), this.backgroundColorInput = O('<input type="color">'), this.aladinDiv = A, A.appendChild(this.mainDiv), this.imgLayers = /* @__PURE__ */ new Map(), this.backgroundColor = this.aladin.getBackgroundColor();
    let I = this;
    this.unselectAllLayers = () => {
      I.aladin.getImageOverlays().forEach((E) => {
        let t = I.imgLayers.get(E).headerDiv[0];
        t.style.backgroundColor = "#f2f2f2";
        let s = t.querySelector(".aladin-layer-header");
        s.style.backgroundColor = "#f2f2f2";
      });
    }, this.selectLayer = (E) => {
      const o = E.layer.layer;
      let t = E.headerDiv[0];
      t.style.backgroundColor = "lightgray";
      let s = t.querySelector(".aladin-layer-header");
      s.style.backgroundColor = "lightgray", I.aladin.setActiveHiPSLayer(o);
    }, this.updateSelectedLayer = () => {
      I.unselectAllLayers();
      const E = I.aladin.getActiveHiPSLayer();
      let o = I.imgLayers.get(E);
      I.selectLayer(o);
    }, this._createComponent(), this._addListeners();
  }
  _onAddCatalogue() {
    if (!this.catalogSelector) {
      let A = function(g, C) {
        if (g == "coneSearch") {
          let I;
          if (C.baseURL.includes("/vizier."))
            I = zA.catalogFromVizieR(
              C.id.replace("CDS/", ""),
              C.ra + " " + C.dec,
              C.radiusDeg,
              { limit: C.limit, onClick: "showTable" }
            );
          else {
            let E = C.baseURL;
            E.endsWith("?") || (E += "?"), E += "RA=" + C.ra + "&DEC=" + C.dec + "&SR=" + C.radiusDeg, I = zA.catalogFromURL(E, { limit: C.limit, onClick: "showTable" });
          }
          this.aladin.addCatalog(I);
        } else if (g == "hips") {
          const I = zA.catalogHiPS(C.hipsURL, { onClick: "showTable", name: C.id });
          this.aladin.addCatalog(I);
        } else if (g == "votable") {
          let I = zA.catalogFromURL(C.url, { onClick: "showTable" });
          this.aladin.addCatalog(I);
        }
      };
      this.catalogSelector = new SE(this.aladinDiv, this.aladin, A);
    }
    this.catalogSelector.show();
  }
  _createComponent() {
    let A = this, g = O(this.mainDiv);
    g.empty(), g.append(
      '<a class="aladin-closeBtn">&times;</a><div class="aladin-box-title">Stack</div>'
    ), g.append('<div class="aladin-box-separator"></div><div class="aladin-label">Image layers</div>'), this.imgLayers.size > 1 && (g.append(
      '<div class="aladin-label" style="font-size: 12px">Overlays</div>'
    ), Array.from(A.aladin.getImageOverlays()).reverse().forEach((T) => {
      let v = A.imgLayers.get(T);
      v && v.layer.layer !== "base" && v.attachTo(g);
    })), g.append(
      '<div class="aladin-label" style="font-size: 12px">Base</div>'
    ), this.imgLayers.has("base") && this.imgLayers.get("base").attachTo(g), g.append(
      '<button class="aladin-btn add-layer-hips" type="button" title="Add a full survey (i.e. a HiPS)">Add survey</button><button class="aladin-btn add-layer-image" type="button" title="Add a single image (only FITS file supported)">Open image 📂</button>'
    ), O(this.mainDiv).find(".add-layer-hips").on("click", function() {
      A.aladin.addNewImageLayer();
    }), O(this.mainDiv).find(".add-layer-image").on("click", function() {
      let T = document.createElement("input");
      T.type = "file", T.onchange = (v) => {
        Array.from(T.files).forEach((a) => {
          const iA = URL.createObjectURL(a), IA = a.name, z = A.aladin.createImageFITS(
            iA,
            IA,
            void 0,
            (oA, hA, yA, rA) => {
              A.aladin.gotoRaDec(oA, hA), A.aladin.setFoV(yA * 1.1);
            },
            void 0
          );
          A.aladin.setOverlayImageLayer(z, IA);
        });
      }, T.click();
    }), g.append('<div class="aladin-box-separator"></div><div class="aladin-label">Overlay layers</div>');
    for (var C = this.aladin.getOverlays(), I = '<ul class="aladin-overlay-list">', E = C.length - 1; E >= 0; E--) {
      var o = C[E], t = o.name, l = "";
      o.isShowing && (l = 'checked="checked"');
      var s = "", M = "";
      if (o.type == "catalog" || o.type == "progressivecat") {
        var N = o.getSources().length;
        s = N + " source" + (N > 1 ? "s" : ""), M = jA.SVG_ICONS.CATALOG;
      } else
        o.type == "moc" ? (s = "Coverage: " + (100 * o.skyFraction()).toFixed(3) + " % of sky", M = jA.SVG_ICONS.MOC) : o.type == "overlay" && (M = jA.SVG_ICONS.OVERLAY);
      var k = O("<div></div>").css("color", o.color).css("color"), U = _A.getLabelColorForBackground(k), Y = window.btoa(M.replace(/FILLCOLOR/g, o.color));
      I += `<li><div class="aladin-stack-icon" style='background-image: url("data:image/svg+xml;base64,` + Y + `");'></div>`, I += '<input class="aladin-input" type="checkbox" ' + l + ' id="aladin_lite_' + o.uuid + '"></input><label for="aladin_lite_' + o.uuid + '" class="aladin-layer-label" style="background: ' + o.color + "; color:" + U + ';" title="' + s + '">' + t + "</label>", I += ' <button class="aladin-btn-small aladin-delete-graphic-layer" type="button" title="Delete this layer" data-uuid="' + o.uuid + '" style="font-size: 10px!important; vertical-align: text-bottom!important; background-color: unset!important;">❌</button>', I += "</li>";
    }
    I += "</ul>", I += '<button class="aladin-btn my-1 catalogue-selector" type="button">Add catalogue</button>', g.append(I), g.find(".aladin-delete-graphic-layer").click(function() {
      const T = A.aladin.findLayerByUUID(O(this).data("uuid"));
      A.aladin.removeLayer(T);
    }), g.find(".catalogue-selector").on("click", () => A._onAddCatalogue()), g.append('<div class="aladin-blank-separator"></div>');
    var l = "";
    A.aladin.isReticleDisplayed() && (l = 'checked="checked"');
    var d = O('<input class="aladin-input" type="checkbox" ' + l + ' id="displayReticle" />');
    g.append(d).append('<label for="displayReticle">Reticle</label><br/>'), d.change(function() {
      A.aladin.showReticle(O(this).is(":checked"));
    }), l = "", A.aladin.isHpxGridDisplayed() && (l = 'checked="checked"');
    var q = O('<input class="aladin-input" type="checkbox" ' + l + ' id="displayHpxGrid"/>');
    g.append(q).append('<label for="displayHpxGrid">HEALPix grid</label><br/>'), q.change(function() {
      A.aladin.showHealpixGrid(O(this).is(":checked"));
    }), g.append('<div class="aladin-box-separator"></div><div class="aladin-label">Background color</div>'), g.append(this.backgroundColorInput), this.backgroundColorInput.on("input", () => {
      A.backgroundColor = this.backgroundColorInput.val(), A.aladin.setBackgroundColor(_A.hexToRgb(A.backgroundColor));
    }), g.append('<div class="aladin-box-separator"></div><div class="aladin-label">Tools</div>');
    var W = O('<button class="aladin-btn" type="button">Export view as PNG</button>');
    g.append(W), W.click(function() {
      A.aladin.exportAsPNG();
    }), g.find(".aladin-closeBtn").click(function() {
      return A.aladin.hideBoxes(), !1;
    }), O(this.mainDiv).find("ul input").change(function() {
      var T = O(this).attr("id").substr(12);
      const v = A.aladin.findLayerByUUID(T);
      O(this).is(":checked") ? v.show() : v.hide();
    });
  }
  _addListeners() {
    let A = this;
    this.aladin.aladinDiv.addEventListener("remove-layer", (g) => {
      const C = g.detail;
      A.aladin.removeImageLayer(C);
    }), this.aladin.aladinDiv.addEventListener("select-layer", (g) => {
      let C = g.detail;
      A.unselectAllLayers(), A.selectLayer(C);
    }), CA.BACKGROUND_COLOR_CHANGED.listenedBy(this.aladin.aladinDiv, function(g) {
      const C = g.detail.color;
      let I = A.mainDiv.querySelector('input[type="color"]'), E = _A.rgbToHex(C.r, C.g, C.b);
      I.value = E, A.backgroundColor = C;
    }), CA.HIPS_LAYER_ADDED.listenedBy(this.aladin.aladinDiv, function(g) {
      const C = g.detail.layer, I = new iQ(A.aladin, C);
      A.imgLayers.set(C.layer, I), A._createComponent(), A.updateSelectedLayer();
    }), CA.HIPS_LAYER_RENAMED.listenedBy(this.aladin.aladinDiv, function(g) {
      const C = g.detail.layer, I = g.detail.newLayer, E = A.imgLayers.get(C);
      A.imgLayers.delete(C), A.imgLayers.set(I, new iQ(A.aladin, E.layer)), A._createComponent(), A.updateSelectedLayer();
    }), CA.HIPS_LAYER_SWAP.listenedBy(this.aladin.aladinDiv, function(g) {
      const C = g.detail.firstLayer, I = g.detail.secondLayer, E = A.imgLayers.get(C), o = A.imgLayers.get(I);
      A.imgLayers.set(I, E), A.imgLayers.set(C, o), A._createComponent(), A.updateSelectedLayer();
    }), CA.HIPS_LAYER_REMOVED.listenedBy(this.aladin.aladinDiv, function(g) {
      const C = g.detail.layer;
      let I = A.imgLayers.get(C);
      I.children ? I.children.forEach((E) => {
        E.destroy(), A.imgLayers.delete(E.layer);
      }) : (I.destroy(), A.imgLayers.delete(C)), A._createComponent(), A.imgLayers.length > 0 && A.updateSelectedLayer();
    }), CA.GRAPHIC_OVERLAY_LAYER_ADDED.listenedBy(this.aladin.aladinDiv, function(g) {
      A._createComponent();
    }), CA.GRAPHIC_OVERLAY_LAYER_REMOVED.listenedBy(this.aladin.aladinDiv, function(g) {
      A._createComponent();
    });
  }
  show() {
    this.mainDiv.style.display = "initial";
  }
  hide() {
    this.mainDiv.style.display = "none";
  }
}
class dE {
  // Constructor
  constructor(A, g, C) {
    this.aladin = g, this.view = C, this.isChecked = !1, this.mainDiv = document.createElement("div"), this.mainDiv.style.display = "none", this.mainDiv.classList.add("aladin-box", "aladin-layerBox", "aladin-cb-list"), this.aladinDiv = A, A.appendChild(this.mainDiv), this._createComponent(), this._addListeners();
  }
  _createComponent() {
    let A = this, g = O(this.mainDiv);
    g.empty(), g.append(
      '<a class="aladin-closeBtn">&times;</a>'
    );
    let C = O('<div class="aladin-label">Coo grid options</div>'), I = O('<div class="layer-options"><table><tbody><tr><td>Color</td><td><input type="color" value="#00ff00"></td></tr><tr><td>Opacity</td><td><input class="aladin-input opacity" value="0.5" type="range" min="0.0" max="1.0" step="0.05"></td></tr><tr><td>Thickness</td><td><input class="aladin-input thickness" value="3.0" type="range" min="0.5" max="10.0" step="0.01"></td></tr><tr><td>Label size</td><td><input class="aladin-input label-size" type="range" value="15" min="5" max="30" step="0.01"></td></tr></table></div>');
    g.append(C).append(I);
    let E = I.find('input[type="color"]'), o = I.find(".opacity"), t = I.find(".thickness");
    t.on("input", () => {
      const N = +t.val();
      A.view.setGridConfig({
        thickness: N
      });
    });
    let s = function() {
      let N = _A.hexToRgb(E.val()), k = o.val();
      A.view.setGridConfig({
        color: { r: N.r / 255, g: N.g / 255, b: N.b / 255 },
        opacity: parseFloat(k)
      });
    };
    E.on("input", s), o.on("input", s);
    let M = I.find(".label-size");
    M.on("input", function() {
      const N = +M.val();
      A.view.setGridConfig({
        labelSize: N
      });
    }), CA.COO_GRID_ENABLED.listenedBy(A.aladinDiv, function() {
      A.isChecked = !A.isChecked;
    }), CA.COO_GRID_DISABLED.listenedBy(A.aladinDiv, function() {
      A.isChecked = !A.isChecked;
    }), CA.COO_GRID_UPDATED.listenedBy(A.aladinDiv, function(N) {
      let k = N.detail.opacity;
      o.val() != k && o.val(k);
      let U = N.detail.color, Y = _A.rgbToHex(Math.round(255 * U.r), Math.round(255 * U.g), Math.round(255 * U.b));
      E.val() != Y && E.val(Y);
    }), g.find(".aladin-closeBtn").click(function() {
      return A.aladin.hideBoxes(), !1;
    });
  }
  _addListeners() {
  }
  show() {
    this.mainDiv.style.display = "block";
  }
  hide() {
    this.mainDiv.style.display = "none";
  }
}
class HE {
  constructor(A) {
    this.aladin = A, this.isShowing = !1;
  }
  _hideMenu(A) {
    this.contextMenuUl.remove(), document.removeEventListener("click", this._hideMenu), window.removeEventListener("resize", this._hideOnResize), this.isShowing = !1;
  }
  _hideOnResize() {
    this._hideMenu(!0);
  }
  _attachOption(A, g, C) {
    const I = document.createElement("li");
    if (I.className = "aladin-context-menu-item", g.label == "Copy position")
      try {
        const o = this.aladin.pix2world(C.x, C.y), t = new $A(o[0], o[1], 6);
        let s;
        this.aladin.view.cooFrame == RA.J2000 ? s = t.format("s/") : (this.aladin.view.cooFrame == RA.J2000d, s = t.format("d/")), I.innerHTML = "<span>" + s + "</span>";
      } catch {
        I.innerHTML = "<span></span>";
      }
    else
      I.innerHTML = "<span>" + g.label + "</span>";
    g.subMenu && g.subMenu.length > 0 && (I.innerHTML += '<span style="position: absolute; right: 4px;">▶</span>');
    const E = this;
    if (I.addEventListener("click", (o) => {
      o.stopPropagation(), (!g.subMenu || g.subMenu.length === 0) && (g.label == "Copy position" ? g.action(o) : g.action(this.event), E._hideMenu(!0));
    }), A.appendChild(I), g.subMenu && g.subMenu.length) {
      const o = document.createElement("ul");
      o.className = "aladin-context-sub-menu", I.appendChild(o), g.subMenu.forEach((t) => this._attachOption(o, t));
    }
  }
  _showMenu(A) {
    this.contextMenuUl.className = "aladin-context-menu", this.contextMenuUl.innerHTML = "";
    const g = Z.relMouseCoords(this.aladin.view.imageCanvas, A);
    this.menuOptions.forEach((M) => this._attachOption(this.contextMenuUl, M, g)), document.body.appendChild(this.contextMenuUl);
    const { innerWidth: C, innerHeight: I } = window, { offsetWidth: E, offsetHeight: o } = this.contextMenuUl;
    let t = 0, s = 0;
    this.event = A, A.clientX >= C / 2 && this.contextMenuUl.classList.add("left"), A.clientY >= I / 2 && this.contextMenuUl.classList.add("top"), A.clientX >= C - E && (t = "-100%"), A.clientY >= I - o && (s = "-100%"), this.contextMenuUl.style.left = A.clientX + "px", this.contextMenuUl.style.top = A.clientY + "px", this.contextMenuUl.style.transform = `translate(${t}, ${s})`, document.addEventListener("click", () => this._hideMenu(!0)), window.addEventListener("resize", this._hideOnResize), this.isShowing = !0;
  }
  attachTo(A, g) {
    this.contextMenuUl = document.createElement("ul"), this.menuOptions = g;
  }
}
let pE = function() {
  function B(g, C, I, E, o = void 0, t = void 0) {
    this.view = I, this.wasm = I.wasm, this.layer = null, this.added = !1, this.subtype = "fits", this.url = g.toString(), this.id = g.toString(), this.name = C, this.imgFormat = "fits", this.properties = {
      formats: ["fits"]
    }, this.successCallback = o, this.errorCallback = t, E && (E.stretch = E.stretch || "asinh"), this.colorCfg = new OI(E);
    let s = this;
    A(s), Fg.update(s), this.query = Promise.resolve(s);
  }
  B.prototype.isReady = function() {
    return this.added;
  }, B.prototype.setOpacity = function(g) {
    let C = this;
    A(C, () => {
      C.colorCfg.setOpacity(g);
    });
  }, B.prototype.setBlendingConfig = function(g = !1) {
    A(this, () => {
      this.colorCfg.setBlendingConfig(g);
    });
  }, B.prototype.setColormap = function(g, C) {
    A(this, () => {
      this.colorCfg.setColormap(g, C);
    });
  }, B.prototype.setCuts = function(g, C) {
    A(this, () => {
      this.colorCfg.setCuts(g, C);
    });
  }, B.prototype.setGamma = function(g) {
    A(this, () => {
      this.colorCfg.setGamma(g);
    });
  }, B.prototype.setSaturation = function(g) {
    A(this, () => {
      this.colorCfg.setSaturation(g);
    });
  }, B.prototype.setBrightness = function(g) {
    A(this, () => {
      this.colorCfg.setBrightness(g);
    });
  }, B.prototype.setContrast = function(g) {
    A(this, () => {
      this.colorCfg.setContrast(g);
    });
  }, B.prototype.metadata = function() {
    return {
      ...this.colorCfg.get(),
      longitudeReversed: !1,
      imgFormat: this.imgFormat
    };
  };
  var A = function(g, C = void 0) {
    C && C();
    try {
      if (g.added) {
        const I = g.metadata();
        g.wasm.setImageMetadata(g.layer, I), CA.HIPS_LAYER_CHANGED.dispatchedTo(g.view.aladinDiv, { layer: g });
      }
    } catch (I) {
      console.error(I);
    }
  };
  return B.prototype.add = function(g) {
    this.layer = g;
    let C = this;
    return C.wasm.addImageFITS({
      layer: C.layer,
      url: C.url,
      meta: C.metadata()
    }).then((E) => {
      C.added = !0, C.children = [];
      let o = 0;
      return E.forEach((t) => {
        let s = new B(
          t.url,
          C.name + "_ext_" + o.toString(),
          C.view,
          null,
          null,
          null
        );
        s.layer = t.layer, s.added = !0, s.colorCfg = Z.deepCopy(C.colorCfg), s.setCuts(t.automatic_min_cut, t.automatic_max_cut), s.ra = t.centered_fov.ra, s.dec = t.centered_fov.dec, s.fov = t.centered_fov.fov, C.ra || (C.ra = s.ra), C.dec || (C.dec = s.dec), C.fov || (C.fov = s.fov), C.children.push(s), o += 1;
      }), C.successCallback && C.successCallback(
        C.children[0].ra,
        C.children[0].dec,
        C.children[0].fov,
        C.children[0]
      ), C;
    }).catch((E) => (window.alert(E + ". See the console for more logging details. It may be possible CORS headers have not been set in the server where you want to download the file. If it is the case, try to manually download the FITS file first and then open it into aladin lite (e.g. by a drag and drop)"), C.errorCallback && C.errorCallback(), C.view.removeImageLayer(g), Promise.reject(E)));
  }, B.prototype.toggle = function() {
    this.colorCfg.getOpacity() != 0 ? this.colorCfg.setOpacity(0) : this.colorCfg.setOpacity(this.prevOpacity);
  }, B.prototype.isPlanetaryBody = function() {
    return !1;
  }, B.prototype.focusOn = function() {
    this.added && (this.view.aladin.gotoRaDec(this.ra, this.dec), this.view.aladin.setFoV(this.fov));
  }, B.prototype.setAlpha = B.prototype.setOpacity, B.prototype.setColorCfg = function(g) {
    A(this, () => {
      this.colorCfg = g;
    });
  }, B.prototype.getColorCfg = function() {
    return this.colorCfg;
  }, B.prototype.getOpacity = function() {
    return this.colorCfg.getOpacity();
  }, B.prototype.getAlpha = B.prototype.getOpacity, B.prototype.readPixel = function(g, C) {
    return this.wasm.readPixel(g, C, this.layer);
  }, B;
}(), fE = function() {
  let B = {};
  return B.getDefaultActions = function(A) {
    return [
      {
        label: "Copy position",
        action(g) {
          var C = document.createRange();
          C.selectNode(g.target), window.getSelection().removeAllRanges(), window.getSelection().addRange(C);
          try {
            let E = document.execCommand("copy") ? "successful" : "unsuccessful";
          } catch {
            console.error("Oops, unable to copy to clipboard");
          }
          window.getSelection().removeAllRanges();
        }
      },
      {
        label: "Take snapshot",
        action(g) {
          A.exportAsPNG();
        }
      },
      {
        label: "Add",
        subMenu: [
          {
            label: "New image layer",
            action(g) {
              A.addNewImageLayer();
            }
          },
          {
            label: "New catalogue layer",
            action(g) {
              A.stack._onAddCatalogue();
            }
          }
        ]
      },
      {
        label: "Load local file",
        subMenu: [
          {
            label: "FITS image",
            action(g) {
              let C = document.createElement("input");
              C.type = "file", C.onchange = (I) => {
                Array.from(C.files).forEach((o) => {
                  const t = URL.createObjectURL(o), s = o.name, M = A.createImageFITS(
                    t,
                    s,
                    void 0,
                    (N, k, U, Y) => {
                      A.gotoRaDec(N, k), A.setFoV(U * 1.1);
                    },
                    void 0
                  );
                  A.setOverlayImageLayer(M, s);
                });
              }, C.click();
            }
          },
          {
            label: "FITS MOC",
            action(g) {
              let C = document.createElement("input");
              C.type = "file", C.onchange = (I) => {
                Array.from(C.files).forEach((o) => {
                  const t = URL.createObjectURL(o);
                  let s = zA.MOCFromURL(t, { name: o.name, fill: !0, opacity: 0.4 });
                  A.addMOC(s);
                });
              }, C.click();
            }
          },
          {
            label: "VOTable",
            action(g) {
              let C = document.createElement("input");
              C.type = "file", C.onchange = (I) => {
                Array.from(C.files).forEach((o) => {
                  const t = URL.createObjectURL(o);
                  let s = zA.catalogFromURL(t, { name: o.name, onClick: "showTable" }, null, !1);
                  A.addCatalog(s);
                });
              }, C.click();
            }
          }
        ]
      },
      {
        label: "What is this?",
        action(g) {
          GQ(A.view, g);
        }
      },
      {
        label: "HiPS2FITS cutout",
        action(g) {
          const C = A;
          let I = "https://alasky.cds.unistra.fr/hips-image-services/hips2fits#", E = C.getRaDec(), o = Math.max.apply(null, C.getFov()), t = C.getBaseImageLayer().id, s = C.getProjectionName();
          I += "ra=" + E[0] + "&dec=" + E[1] + "&fov=" + o + "&projection=" + s + "&hips=" + encodeURIComponent(t), window.open(I, "_blank");
        }
      },
      {
        label: "Select sources",
        action(g) {
          A.select();
        }
      }
    ];
  }, B;
}(), VA = function() {
  var B = function(I, E) {
    if (O(I).length == 0)
      return;
    this.wasm = null;
    const o = this;
    if (E === void 0 && (E = this.getOptionsFromQueryString()), E = E || {}, "zoom" in E) {
      var t = E.zoom;
      delete E.zoom, E.fov = t;
    }
    var s = {};
    for (var M in B.DEFAULT_OPTIONS)
      E[M] !== void 0 ? s[M] = E[M] : s[M] = B.DEFAULT_OPTIONS[M];
    for (var M in E)
      B.DEFAULT_OPTIONS[M] === void 0 && (s[M] = E[M]);
    this.options = s, O("<style type='text/css'> .aladin-reticleColor { color: " + this.options.reticleColor + "; font-weight:bold;} </style>").appendTo(I), this.aladinDiv = I, this.reduceDeformations = !0, O(I).addClass("aladin-container");
    let N = RA.fromString(s.cooFrame, RA.J2000);
    const k = O('<div class="aladin-location">' + (s.showFrame ? '<select class="aladin-selector aladin-frameChoice"><option value="' + RA.J2000.label + '" ' + (N == RA.J2000 ? 'selected="selected"' : "") + '>J2000</option><option value="' + RA.J2000d.label + '" ' + (N == RA.J2000d ? 'selected="selected"' : "") + '>J2000d</option><option value="' + RA.GAL.label + '" ' + (N == RA.GAL ? 'selected="selected"' : "") + ">GAL</option></select>" : "") + '<span class="aladin-clipboard" title="Copy coordinates to clipboard"></span><span class="aladin-location-text"></span></div>').appendTo(I), U = k.find(".aladin-clipboard");
    U.hide(), U.click(function() {
      o.copyCoordinatesToClipboard();
    }), k.mouseenter(function() {
      U.show();
    }), k.mouseleave(function() {
      U.hide();
    });
    var Y = O('<div class="aladin-fov"></div>').appendTo(I);
    s.showZoomControl && O('<div class="aladin-zoomControl"><a href="#" class="zoomPlus" title="Zoom in">+</a><a href="#" class="zoomMinus" title="Zoom out">&ndash;</a></div>').appendTo(I), s.showFullscreenControl && O('<div class="aladin-fullscreenControl aladin-maximize" title="Full screen"></div>').appendTo(I), this.fullScreenBtn = O(I).find(".aladin-fullscreenControl"), this.fullScreenBtn.click(function() {
      o.toggleFullscreen(o.options.realFullscreen);
    }), O(document).on("fullscreenchange webkitfullscreenchange mozfullscreenchange MSFullscreenChange", function(rA) {
      var gA = document.fullscreenElement || document.webkitFullscreenElement || document.mozFullScreenElement || document.msFullscreenElement;
      if (gA == null) {
        o.fullScreenBtn.removeClass("aladin-restore"), o.fullScreenBtn.addClass("aladin-maximize"), o.fullScreenBtn.attr("title", "Full screen"), O(o.aladinDiv).removeClass("aladin-fullscreen");
        var nA = o.callbacksByEventName.fullScreenToggled, og = o.fullScreenBtn.hasClass("aladin-restore");
        typeof nA == "function" && nA(og);
      }
    }), new lE(I), this.boxes = [], this.measurementTable = new RE(I);
    var K = new KE(k.find(".aladin-location-text"));
    this.view = new eI(this, K, Y, N, s.fov), this.cacheSurveys = /* @__PURE__ */ new Map(), this.stack = new qE(this.aladinDiv, this, this.view), this.coogrid = new dE(this.aladinDiv, this, this.view), s.backgroundColor && (this.backgroundColor = s.backgroundColor, this.setBackgroundColor(this.backgroundColor)), this.boxes.push(this.stack), this.boxes.push(this.coogrid);
    var l, d;
    s.gridOptions ? (l = s.gridOptions.color && _A.hexToRgb(s.gridOptions.color), d = s.gridOptions.opacity) : (l = { r: 0, g: 1, b: 0 }, d = 0.5), this.view.setGridConfig({
      color: l,
      opacity: d
    }), s && s.showCooGrid && this.showCooGrid(), s && (s.showProjectionControl === void 0 || s.showProjectionControl === !0) && new JE(I, this);
    let q = s && s.projection || "SIN";
    this.setProjection(q);
    let W = 30;
    if (CA.LOADING_STATE.listenedBy(I, function(rA) {
      let gA = I.querySelector(".aladin-layersControl");
      gA && (rA.detail.loading ? gA.style.backgroundImage = "url(data:image/svg+xml;base64,PD94bWwgdmVyc2lvbj0iMS4wIiBzdGFuZGFsb25lPSJubyI/Pgo8IURPQ1RZUEUgc3ZnIFBVQkxJQyAiLS8vVzNDLy9EVEQgU1ZHIDEuMS8vRU4iICJodHRwOi8vd3d3LnczLm9yZy9HcmFwaGljcy9TVkcvMS4xL0RURC9zdmcxMS5kdGQiPgo8c3ZnIHdpZHRoPSI0MHB4IiBoZWlnaHQ9IjQwcHgiIHZpZXdCb3g9IjAgMCA0MCA0MCIgdmVyc2lvbj0iMS4xIiB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHhtbG5zOnhsaW5rPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5L3hsaW5rIiB4bWw6c3BhY2U9InByZXNlcnZlIiBzdHlsZT0iZmlsbC1ydWxlOmV2ZW5vZGQ7Y2xpcC1ydWxlOmV2ZW5vZGQ7c3Ryb2tlLWxpbmVqb2luOnJvdW5kO3N0cm9rZS1taXRlcmxpbWl0OjEuNDE0MjE7IiB4PSIwcHgiIHk9IjBweCI+CiAgICA8ZGVmcz4KICAgICAgICA8c3R5bGUgdHlwZT0idGV4dC9jc3MiPjwhW0NEQVRBWwogICAgICAgICAgICBALXdlYmtpdC1rZXlmcmFtZXMgc3BpbiB7CiAgICAgICAgICAgICAgZnJvbSB7CiAgICAgICAgICAgICAgICAtd2Via2l0LXRyYW5zZm9ybTogcm90YXRlKDBkZWcpCiAgICAgICAgICAgICAgfQogICAgICAgICAgICAgIHRvIHsKICAgICAgICAgICAgICAgIC13ZWJraXQtdHJhbnNmb3JtOiByb3RhdGUoLTM1OWRlZykKICAgICAgICAgICAgICB9CiAgICAgICAgICAgIH0KICAgICAgICAgICAgQGtleWZyYW1lcyBzcGluIHsKICAgICAgICAgICAgICBmcm9tIHsKICAgICAgICAgICAgICAgIHRyYW5zZm9ybTogcm90YXRlKDBkZWcpCiAgICAgICAgICAgICAgfQogICAgICAgICAgICAgIHRvIHsKICAgICAgICAgICAgICAgIHRyYW5zZm9ybTogcm90YXRlKC0zNTlkZWcpCiAgICAgICAgICAgICAgfQogICAgICAgICAgICB9CiAgICAgICAgICAgIHN2ZyB7CiAgICAgICAgICAgICAgICAtd2Via2l0LXRyYW5zZm9ybS1vcmlnaW46IDUwJSA1MCU7CiAgICAgICAgICAgICAgICAtd2Via2l0LWFuaW1hdGlvbjogc3BpbiAxLjVzIGxpbmVhciBpbmZpbml0ZTsKICAgICAgICAgICAgICAgIC13ZWJraXQtYmFja2ZhY2UtdmlzaWJpbGl0eTogaGlkZGVuOwogICAgICAgICAgICAgICAgYW5pbWF0aW9uOiBzcGluIDEuNXMgbGluZWFyIGluZmluaXRlOwogICAgICAgICAgICB9CiAgICAgICAgXV0+PC9zdHlsZT4KICAgIDwvZGVmcz4KICAgIDxnIGlkPSJvdXRlciI+CiAgICAgICAgPGc+CiAgICAgICAgICAgIDxwYXRoIGQ9Ik0yMCwwQzIyLjIwNTgsMCAyMy45OTM5LDEuNzg4MTMgMjMuOTkzOSwzLjk5MzlDMjMuOTkzOSw2LjE5OTY4IDIyLjIwNTgsNy45ODc4MSAyMCw3Ljk4NzgxQzE3Ljc5NDIsNy45ODc4MSAxNi4wMDYxLDYuMTk5NjggMTYuMDA2MSwzLjk5MzlDMTYuMDA2MSwxLjc4ODEzIDE3Ljc5NDIsMCAyMCwwWiIgc3R5bGU9ImZpbGw6YmxhY2s7Ii8+CiAgICAgICAgPC9nPgogICAgICAgIDxnPgogICAgICAgICAgICA8cGF0aCBkPSJNNS44NTc4Niw1Ljg1Nzg2QzcuNDE3NTgsNC4yOTgxNSA5Ljk0NjM4LDQuMjk4MTUgMTEuNTA2MSw1Ljg1Nzg2QzEzLjA2NTgsNy40MTc1OCAxMy4wNjU4LDkuOTQ2MzggMTEuNTA2MSwxMS41MDYxQzkuOTQ2MzgsMTMuMDY1OCA3LjQxNzU4LDEzLjA2NTggNS44NTc4NiwxMS41MDYxQzQuMjk4MTUsOS45NDYzOCA0LjI5ODE1LDcuNDE3NTggNS44NTc4Niw1Ljg1Nzg2WiIgc3R5bGU9ImZpbGw6cmdiKDIxMCwyMTAsMjEwKTsiLz4KICAgICAgICA8L2c+CiAgICAgICAgPGc+CiAgICAgICAgICAgIDxwYXRoIGQ9Ik0yMCwzMi4wMTIyQzIyLjIwNTgsMzIuMDEyMiAyMy45OTM5LDMzLjgwMDMgMjMuOTkzOSwzNi4wMDYxQzIzLjk5MzksMzguMjExOSAyMi4yMDU4LDQwIDIwLDQwQzE3Ljc5NDIsNDAgMTYuMDA2MSwzOC4yMTE5IDE2LjAwNjEsMzYuMDA2MUMxNi4wMDYxLDMzLjgwMDMgMTcuNzk0MiwzMi4wMTIyIDIwLDMyLjAxMjJaIiBzdHlsZT0iZmlsbDpyZ2IoMTMwLDEzMCwxMzApOyIvPgogICAgICAgIDwvZz4KICAgICAgICA8Zz4KICAgICAgICAgICAgPHBhdGggZD0iTTI4LjQ5MzksMjguNDkzOUMzMC4wNTM2LDI2LjkzNDIgMzIuNTgyNCwyNi45MzQyIDM0LjE0MjEsMjguNDkzOUMzNS43MDE5LDMwLjA1MzYgMzUuNzAxOSwzMi41ODI0IDM0LjE0MjEsMzQuMTQyMUMzMi41ODI0LDM1LjcwMTkgMzAuMDUzNiwzNS43MDE5IDI4LjQ5MzksMzQuMTQyMUMyNi45MzQyLDMyLjU4MjQgMjYuOTM0MiwzMC4wNTM2IDI4LjQ5MzksMjguNDkzOVoiIHN0eWxlPSJmaWxsOnJnYigxMDEsMTAxLDEwMSk7Ii8+CiAgICAgICAgPC9nPgogICAgICAgIDxnPgogICAgICAgICAgICA8cGF0aCBkPSJNMy45OTM5LDE2LjAwNjFDNi4xOTk2OCwxNi4wMDYxIDcuOTg3ODEsMTcuNzk0MiA3Ljk4NzgxLDIwQzcuOTg3ODEsMjIuMjA1OCA2LjE5OTY4LDIzLjk5MzkgMy45OTM5LDIzLjk5MzlDMS43ODgxMywyMy45OTM5IDAsMjIuMjA1OCAwLDIwQzAsMTcuNzk0MiAxLjc4ODEzLDE2LjAwNjEgMy45OTM5LDE2LjAwNjFaIiBzdHlsZT0iZmlsbDpyZ2IoMTg3LDE4NywxODcpOyIvPgogICAgICAgIDwvZz4KICAgICAgICA8Zz4KICAgICAgICAgICAgPHBhdGggZD0iTTUuODU3ODYsMjguNDkzOUM3LjQxNzU4LDI2LjkzNDIgOS45NDYzOCwyNi45MzQyIDExLjUwNjEsMjguNDkzOUMxMy4wNjU4LDMwLjA1MzYgMTMuMDY1OCwzMi41ODI0IDExLjUwNjEsMzQuMTQyMUM5Ljk0NjM4LDM1LjcwMTkgNy40MTc1OCwzNS43MDE5IDUuODU3ODYsMzQuMTQyMUM0LjI5ODE1LDMyLjU4MjQgNC4yOTgxNSwzMC4wNTM2IDUuODU3ODYsMjguNDkzOVoiIHN0eWxlPSJmaWxsOnJnYigxNjQsMTY0LDE2NCk7Ii8+CiAgICAgICAgPC9nPgogICAgICAgIDxnPgogICAgICAgICAgICA8cGF0aCBkPSJNMzYuMDA2MSwxNi4wMDYxQzM4LjIxMTksMTYuMDA2MSA0MCwxNy43OTQyIDQwLDIwQzQwLDIyLjIwNTggMzguMjExOSwyMy45OTM5IDM2LjAwNjEsMjMuOTkzOUMzMy44MDAzLDIzLjk5MzkgMzIuMDEyMiwyMi4yMDU4IDMyLjAxMjIsMjBDMzIuMDEyMiwxNy43OTQyIDMzLjgwMDMsMTYuMDA2MSAzNi4wMDYxLDE2LjAwNjFaIiBzdHlsZT0iZmlsbDpyZ2IoNzQsNzQsNzQpOyIvPgogICAgICAgIDwvZz4KICAgICAgICA8Zz4KICAgICAgICAgICAgPHBhdGggZD0iTTI4LjQ5MzksNS44NTc4NkMzMC4wNTM2LDQuMjk4MTUgMzIuNTgyNCw0LjI5ODE1IDM0LjE0MjEsNS44NTc4NkMzNS43MDE5LDcuNDE3NTggMzUuNzAxOSw5Ljk0NjM4IDM0LjE0MjEsMTEuNTA2MUMzMi41ODI0LDEzLjA2NTggMzAuMDUzNiwxMy4wNjU4IDI4LjQ5MzksMTEuNTA2MUMyNi45MzQyLDkuOTQ2MzggMjYuOTM0Miw3LjQxNzU4IDI4LjQ5MzksNS44NTc4NloiIHN0eWxlPSJmaWxsOnJnYig1MCw1MCw1MCk7Ii8+CiAgICAgICAgPC9nPgogICAgPC9nPgo8L3N2Zz4K)" : gA.style.backgroundImage = 'url("data:image/gif;base64,R0lGODlhGQAcAMIAAAAAADQ0NKahocvFxf///wAAAAAAAAAAACH5BAEKAAcALAAAAAAZABwAAANneLoH/hCwyaJ1dDrCuydY1gBfyYUaaZqosq0r+sKxNNP1pe98Hy2OgXBILLZGxWRSBlA6iZjgczrwWa9WIEDA7Xq/R8d3PGaSz97oFs0WYN9wiDZAr9vvYcB9v2fy/3ZqgIN0cYZYCQA7")');
    }), s.showLayersControl) {
      var T = O('<div class="aladin-layersControl-container" style="top: ' + W + 'px" title="Manage layers"><div class="aladin-layersControl"></div></div>');
      T.appendTo(I), s.expandLayersControl && (o.hideBoxes(), o.showLayerBox()), T.click(function() {
        return o.hideBoxes(), o.showLayerBox(), !1;
      }), W += 38;
    }
    if (s.showGotoControl) {
      var T = O('<div class="aladin-gotoControl-container" style="top: ' + W + 'px" title="Go to position"><div class="aladin-gotoControl"></div></div>');
      T.appendTo(I);
      var v = O('<div class="aladin-box aladin-gotoBox"><a class="aladin-closeBtn" style="display: inline-block">&times;</a><form class="aladin-target-form" style="display: inline-block">Go to: <input class="aladin-input" type="text" placeholder="Object name/position" /></form></div>');
      v.appendTo(I), this.boxes.push(v);
      var _ = v.find(".aladin-target-form input");
      _.on("paste keydown", function() {
        O(this).removeClass("aladin-unknownObject");
      }), _.on("change", function() {
        _.blur();
      }), T.click(function() {
        return o.hideBoxes(), _.val(""), _.removeClass("aladin-unknownObject"), v.show(), _.blur(), !1;
      }), v.find(".aladin-closeBtn").click(function() {
        return o.hideBoxes(), _.blur(), !1;
      }), W += 38;
    }
    if (s.showSimbadPointerControl) {
      var T = O('<div class="aladin-simbadPointerControl-container" style="top: ' + W + 'px" title="What is this? Click on an object to identify it."><div class="aladin-simbadPointerControl"></div></div>');
      T.appendTo(I), T.click(function() {
        o.view.setMode(eI.TOOL_SIMBAD_POINTER);
      }), W += 38;
    }
    if (s.showCooGridControl) {
      var T = O('<div class="aladin-cooGridControl-container" style="top: ' + W + 'px" title="Coo grid. Keep the mouse down to see the option panel"><div class="aladin-cooGridControl"></div></div>');
      T.appendTo(I);
      let gA = !1, nA, og, PA = !1;
      T.on("mousedown", function() {
        gA = !0, nA = /* @__PURE__ */ new Date(), og = setInterval(() => {
          if (gA && /* @__PURE__ */ new Date() - nA > 500)
            return PA = !0, o.hideBoxes(), o.showCooGridBox(), !1;
        }, 50);
      }), T.on("mouseup", function() {
        gA = !1, clearInterval(og), /* @__PURE__ */ new Date() - nA < 500 && (PA = !1);
      }), T.on("click", function() {
        if (!PA) {
          let gg = T[0];
          o.cooGridEnabled ? (o.hideCooGrid(), gg.style.background = "rgba(250, 250, 250, 0.8)") : (o.showCooGrid(), gg.style.background = "rgba(250, 250, 250, 1)");
        }
      }), W += 38;
    }
    if (s.showShareControl) {
      var T = O('<div class="aladin-shareControl-container" title="Get link for current view"><div class="aladin-shareControl"></div></div>');
      T.appendTo(I);
      var a = O('<div class="aladin-box aladin-shareBox"><a class="aladin-closeBtn">&times;</a><div style="clear: both;"></div>Link to previewer: <span class="info"></span><input type="text" class="aladin-input aladin-shareInput" /></div>');
      a.appendTo(I), this.boxes.push(a), T.click(function() {
        o.hideBoxes(), a.show();
        var gA = o.getShareURL();
        return a.find(".aladin-shareInput").val(gA).select(), document.execCommand("copy"), !1;
      }), a.find(".aladin-closeBtn").click(function() {
        return o.hideBoxes(), !1;
      }), W += 38;
    }
    if (this.gotoObject(s.target, void 0, { forceAnimation: !1 }), s.log) {
      var iA = E;
      iA.version = B.VERSION, cB.log("startup", iA);
    }
    if (this.showReticle(s.showReticle), s.catalogUrls)
      for (var IA = 0, z = s.catalogUrls.length; IA < z; IA++)
        this.createCatalogFromVOTable(s.catalogUrls[IA]);
    if (s.survey)
      if (Array.isArray(s.survey)) {
        let rA = 0;
        s.survey.forEach((gA) => {
          rA == 0 ? this.setBaseImageLayer(gA) : this.setOverlayImageLayer(gA, Z.uuidv4()), rA++;
        });
      } else
        this.setBaseImageLayer(s.survey);
    else {
      const rA = Math.round(Math.random()), gA = B.DEFAULT_OPTIONS.surveyUrl[rA];
      this.setBaseImageLayer(gA);
    }
    this.view.showCatalog(s.showCatalog);
    var oA = this;
    O(I).find(".aladin-frameChoice").change(function() {
      oA.setFrame(O(this).val());
    }), O(I).find(".aladin-target-form").submit(function() {
      return oA.gotoObject(O(this).find("input").val(), function() {
        O(I).find(".aladin-target-form input").addClass("aladin-unknownObject");
      }), !1;
    });
    var hA = O(I).find(".zoomPlus");
    hA.click(function() {
      return oA.increaseZoom(), !1;
    }), hA.bind("mousedown", function(rA) {
      rA.preventDefault();
    });
    var yA = O(I).find(".zoomMinus");
    yA.click(function() {
      return oA.decreaseZoom(), !1;
    }), yA.bind("mousedown", function(rA) {
      rA.preventDefault();
    }), this.callbacksByEventName = {}, this.view.redraw(), s.fullScreen && o.toggleFullscreen(o.options.realFullscreen), s.showContextMenu && (this.contextMenu = new HE(this), this.contextMenu.attachTo(this.view.catalogCanvas, fE.getDefaultActions(this)));
  };
  B.VERSION = "3.0-beta0", B.JSONP_PROXY = "https://alaskybis.cds.unistra.fr/cgi/JSONProxy", B.URL_PREVIEWER = "https://aladin.cds.unistra.fr/AladinLite/", B.wasmLibs = {}, B.DEFAULT_OPTIONS = {
    surveyUrl: ["https://alaskybis.u-strasbg.fr/DSS/DSSColor", "https://alasky.u-strasbg.fr/DSS/DSSColor"],
    survey: "CDS/P/DSS2/color",
    target: "0 +0",
    cooFrame: "J2000",
    fov: 60,
    backgroundColor: "rgb(0, 0, 0)",
    showReticle: !0,
    showZoomControl: !0,
    showFullscreenControl: !0,
    showLayersControl: !0,
    showGotoControl: !0,
    showSimbadPointerControl: !1,
    showShareControl: !1,
    showContextMenu: !1,
    showCatalog: !0,
    // TODO: still used ??
    showFrame: !0,
    fullScreen: !1,
    reticleColor: "rgb(178, 50, 178)",
    reticleSize: 22,
    log: !0,
    allowFullZoomout: !1,
    realFullscreen: !1,
    showAllskyRing: !1,
    allskyRingColor: "#c8c8ff",
    allskyRingWidth: 8,
    pixelateCanvas: !0
  }, B.prototype.copyCoordinatesToClipboard = function() {
    let I = this.view.location.$div[0];
    var E = document.createRange();
    E.selectNode(I), window.getSelection().removeAllRanges(), window.getSelection().addRange(E);
    try {
      let t = document.execCommand("copy") ? "successful" : "unsuccessful";
      console.log("Copying text command was " + t);
    } catch {
      console.log("Oops, unable to copy");
    }
    window.getSelection().removeAllRanges();
  }, B.prototype.toggleFullscreen = function(I) {
    let E = this;
    I = !!I, this.fullScreenBtn.toggleClass("aladin-maximize aladin-restore");
    var o = this.fullScreenBtn.hasClass("aladin-restore");
    if (this.fullScreenBtn.attr("title", o ? "Restore original size" : "Full screen"), this.aladinDiv.classList.contains("aladin-fullscreen") ? this.aladinDiv.classList.remove("aladin-fullscreen") : this.aladinDiv.classList.add("aladin-fullscreen"), I)
      if (o) {
        var t = this.aladinDiv;
        t.requestFullscreen ? t.requestFullscreen() : t.webkitRequestFullscreen ? t.webkitRequestFullscreen() : t.mozRequestFullScreen ? t.mozRequestFullScreen() : t.msRequestFullscreen && t.msRequestFullscreen();
      } else
        document.exitFullscreen ? document.exitFullscreen() : document.webkitExitFullscreen ? document.webkitExitFullscreen() : document.mozCancelFullScreen ? document.mozCancelFullScreen() : document.webkitExitFullscreen && document.webkitExitFullscreen();
    var s = E.callbacksByEventName.zoomChanged;
    typeof s == "function" && s(E.view.fov);
    var M = E.callbacksByEventName.fullScreenToggled;
    typeof M == "function" && M(o);
  }, B.prototype.getOptionsFromQueryString = function() {
    var I = {}, E = Z.urlParam("target");
    E && (I.target = E);
    var o = Z.urlParam("frame");
    o && RA[o] && (I.frame = o);
    var t = Z.urlParam("survey");
    t && EQ.getSurveyInfoFromId(t) && (I.survey = t);
    var s = Z.urlParam("zoom");
    s && s > 0 && s < 180 && (I.zoom = s);
    var M = Z.urlParam("showReticle");
    M && (I.showReticle = M.toLowerCase() == "true");
    var N = Z.urlParam("cooFrame");
    N && (I.cooFrame = N);
    var k = Z.urlParam("fullScreen");
    return k !== void 0 && (I.fullScreen = k), I;
  }, B.prototype.setFoV = B.prototype.setFov = function(I) {
    this.view.setZoom(I);
  }, B.prototype.adjustFovForObject = function(I) {
    var E = this;
    this.getFovForObject(I, function(o) {
      E.setFoV(o);
    });
  }, B.prototype.getFovForObject = function(I, E) {
    var o = "SELECT galdim_majaxis, V FROM basic JOIN ident ON oid=ident.oidref JOIN allfluxes ON oid=allfluxes.oidref WHERE id='" + I + "'", t = "//simbad.u-strasbg.fr/simbad/sim-tap/sync?query=" + encodeURIComponent(o) + "&request=doQuery&lang=adql&format=json&phase=run", s = Z.getAjaxObject(t, "GET", "json", !1);
    s.done(function(M) {
      var N = 0.06666666666666667, k = N;
      if ("data" in M && M.data.length > 0) {
        var U = Z.isNumber(M.data[0][0]) ? M.data[0][0] / 60 : null, Y = Z.isNumber(M.data[0][1]) ? M.data[0][1] : null;
        U !== null ? k = 2 * U : Y !== null && Y < 10 && (k = 2 * Math.pow(2, 6 - Y / 2) / 60);
      }
      typeof E == "function" && E(k);
    });
  }, B.prototype.setFrame = function(I) {
    if (I) {
      var E = RA.fromString(I, RA.J2000);
      if (E != this.view.cooFrame) {
        this.view.changeFrame(E);
        var o = this.view.aladin.callbacksByEventName.cooFrameChanged;
        typeof o == "function" && o(E.label), O(this.aladinDiv).find(".aladin-frameChoice").val(E.label);
      }
    }
  }, B.prototype.setProjection = function(I) {
    I && (this.view.setProjection(I), CA.PROJECTION_CHANGED.dispatchedTo(this.aladinDiv, { projection: I }));
  }, B.prototype.getProjectionName = function() {
    const I = this;
    let E;
    for (let o in uA)
      if (uA[o].id == I.view.projection.id) {
        E = o;
        break;
      }
    return E;
  }, B.prototype.gotoObject = function(I, E, o) {
    let t, s;
    typeof E == "object" ? (E.hasOwnProperty("success") && (t = E.success), E.hasOwnProperty("error") && (s = E.error)) : typeof E == "function" && (s = E);
    var M = /[a-zA-Z]/.test(I);
    if (M) {
      var k = this;
      (async () => {
        let U;
        if (this.getBaseImageLayer() && (U = await this.getBaseImageLayer().query), this.getBaseImageLayer() === void 0 || !U.isPlanetaryBody())
          nE.resolve(
            I,
            function(Y) {
              const K = Y.Target.Resolver;
              k.view.pointTo(K.jradeg, K.jdedeg, o), typeof t == "function" && t(k.getRaDec());
            },
            function(Y) {
              console && (console.log("Could not resolve object name " + I), console.log(Y)), typeof s == "function" && s();
            }
          );
        else {
          const Y = U.properties.hipsBody;
          FE.resolve(
            I,
            Y,
            function(K) {
              k.view.pointTo(K.lon, K.lat, o), typeof t == "function" && t(k.getRaDec());
            },
            function(K) {
              console && (console.log("Could not resolve object name " + I), console.log(K)), typeof s == "function" && s();
            }
          );
        }
      })();
    } else {
      var N = new $A();
      N.parse(I);
      const [U, Y] = this.wasm.viewToICRSCooSys(N.lon, N.lat);
      this.view.pointTo(U, Y, o), typeof t == "function" && t(this.getRaDec());
    }
  }, B.prototype.gotoPosition = function(I, E) {
    var o;
    this.view.cooFrame == RA.GAL ? o = YE.GalacticToJ2000([I, E]) : o = [I, E], this.view.pointTo(o[0], o[1]);
  };
  var A = function(I) {
    var E = I.animationParams;
    if (!(E == null || !E.running)) {
      var o = (/* @__PURE__ */ new Date()).getTime();
      if (o > E.end) {
        I.gotoRaDec(E.raEnd, E.decEnd), E.complete && E.complete();
        return;
      }
      var t = (o - E.start) / (E.end - E.start), s = C(E.raStart, E.decStart, E.raEnd, E.decEnd, t), M = s[0], N = s[1];
      I.gotoRaDec(M, N), setTimeout(function() {
        A(I);
      }, 10);
    }
  };
  B.prototype.stopAnimation = function() {
    this.zoomAnimationParams && (this.zoomAnimationParams.running = !1), this.animationParams && (this.animationParams.running = !1);
  }, B.prototype.animateToRaDec = function(I, E, o, t) {
    o = o || 5, this.animationParams = null;
    var s = {};
    s.start = (/* @__PURE__ */ new Date()).getTime(), s.end = (/* @__PURE__ */ new Date()).getTime() + 1e3 * o;
    var M = this.getRaDec();
    s.raStart = M[0], s.decStart = M[1], s.raEnd = I, s.decEnd = E, s.complete = t, s.running = !0, this.animationParams = s, A(this);
  };
  var g = function(I) {
    var E = I.zoomAnimationParams;
    if (!(E == null || !E.running)) {
      var o = (/* @__PURE__ */ new Date()).getTime();
      if (o > E.end) {
        I.setFoV(E.fovEnd), E.complete && E.complete();
        return;
      }
      var t = (o - E.start) / (E.end - E.start), s = E.fovStart + (E.fovEnd - E.fovStart) * Math.sqrt(t);
      I.setFoV(s), setTimeout(function() {
        g(I);
      }, 50);
    }
  };
  B.prototype.zoomToFoV = function(I, E, o) {
    E = E || 5, this.zoomAnimationParams = null;
    var t = {};
    t.start = (/* @__PURE__ */ new Date()).getTime(), t.end = (/* @__PURE__ */ new Date()).getTime() + 1e3 * E;
    var s = this.getFov();
    t.fovStart = Math.max(s[0], s[1]), t.fovEnd = I, t.complete = o, t.running = !0, this.zoomAnimationParams = t, g(this);
  };
  function C(U, k, K, Y, s) {
    function M(iA) {
      return iA * Math.PI / 180;
    }
    function N(iA) {
      return iA * 180 / Math.PI;
    }
    var k = M(k), U = M(U), Y = M(Y), K = M(K), l = 2 * Math.asin(
      Math.sqrt(Math.pow(
        Math.sin((k - Y) / 2),
        2
      ) + Math.cos(k) * Math.cos(Y) * Math.pow(Math.sin((U - K) / 2), 2))
    ), d = Math.sin((1 - s) * l) / Math.sin(l), q = Math.sin(s * l) / Math.sin(l), W = d * Math.cos(k) * Math.cos(U) + q * Math.cos(Y) * Math.cos(K), T = d * Math.cos(k) * Math.sin(U) + q * Math.cos(Y) * Math.sin(K), v = d * Math.sin(k) + q * Math.sin(Y), _ = Math.atan2(T, W), a = Math.atan2(v, Math.sqrt(Math.pow(W, 2) + Math.pow(T, 2)));
    return [N(_), N(a)];
  }
  return B.prototype.getRaDec = function() {
    let I = this.wasm.getCenter();
    const E = this.wasm.viewToICRSCooSys(I[0], I[1]);
    return E[0] < 0 ? [E[0] + 360, E[1]] : E;
  }, B.prototype.gotoRaDec = function(I, E) {
    this.view.pointTo(I, E);
  }, B.prototype.showHealpixGrid = function(I) {
    this.view.showHealpixGrid(I);
  }, B.prototype.showSurvey = function(I) {
    this.view.showSurvey(I);
  }, B.prototype.showCatalog = function(I) {
    this.view.showCatalog(I);
  }, B.prototype.showReticle = function(I) {
    this.view.showReticle(I), O("#displayReticle").attr("checked", I);
  }, B.prototype.removeLayers = function() {
    this.view.removeLayers();
  }, B.prototype.addCatalog = function(I) {
    this.view.addCatalog(I), CA.GRAPHIC_OVERLAY_LAYER_ADDED.dispatchedTo(this.aladinDiv, { layer: I });
  }, B.prototype.addOverlay = function(I) {
    this.view.addOverlay(I), CA.GRAPHIC_OVERLAY_LAYER_ADDED.dispatchedTo(this.aladinDiv, { layer: I });
  }, B.prototype.addMOC = function(I) {
    this.view.addMOC(I);
  }, B.prototype.findLayerByUUID = function(I) {
    const E = this.view.allOverlayLayers.filter((o) => o.uuid === I);
    return E.length == 0 ? null : E[0];
  }, B.prototype.removeLayer = function(I) {
    this.view.removeLayer(I);
  }, B.prototype.createImageSurvey = function(I, E, o, t, s, M = {}) {
    let N = this.cacheSurveys.get(I);
    return N ? N = Z.clone(N) : (t && (M.cooFrame = t), s && (M.maxOrder = s), N = { id: I, name: E, rootUrl: o, options: M }, this.cacheSurveys.set(I, N)), new EQ(N.id, N.name, N.rootUrl, this.view, N.options);
  }, B.prototype.createImageFITS = function(I, E, o = {}, t = void 0, s = void 0) {
    try {
      I = new URL(I);
    } catch {
      I = Z.getAbsoluteURL(I), I = new URL(I);
    }
    let M = this.cacheSurveys.get(I);
    return M ? M = Z.clone(M) : (M = { url: I, name: E, options: o, successCallback: t, errorCallback: s }, this.cacheSurveys.set(I, M)), new pE(M.url, M.name, this.view, M.options, M.successCallback, M.errorCallback);
  }, B.prototype.newImageSurvey = function(I, E) {
    const o = I, t = o;
    return this.createImageSurvey(o, t, o, null, null, E);
  }, B.prototype.addNewImageLayer = function() {
    let I = Z.uuidv4();
    this.setOverlayImageLayer("CDS/P/DSS2/color", I);
  }, B.prototype.setImageLayer = function(I) {
    this.setBaseImageLayer(I);
  }, B.prototype.setImageSurvey = B.prototype.setImageLayer, B.prototype.setBackgroundColor = function(I) {
    let E;
    if (typeof I == "string") {
      var I = I.match(/^rgb\((\d+),\s*(\d+),\s*(\d+)\)$/), o = parseInt(I[1]), t = parseInt(I[2]), s = parseInt(I[3]);
      E = { r: o, g: t, b: s };
    } else
      E = I;
    this.backgroundColor = E, CA.AL_USE_WASM.dispatchedTo(document.body, { callback: (M) => {
      M.setBackgroundColor(this.backgroundColor), CA.BACKGROUND_COLOR_CHANGED.dispatchedTo(this.aladinDiv, { color: this.backgroundColor });
    } });
  }, B.prototype.getBackgroundColor = function() {
    return this.backgroundColor;
  }, B.prototype.removeImageLayer = function(I) {
    this.view.removeImageLayer(I);
  }, B.prototype.setBaseImageLayer = function(I) {
    return this.setOverlayImageLayer(I, "base");
  }, B.prototype.getBaseImageLayer = function() {
    return this.view.getImageLayer("base");
  }, B.prototype.setOverlayImageLayer = function(I, E = "overlay") {
    let o;
    if (typeof I == "string") {
      const t = I, s = t;
      o = this.createImageSurvey(t, s, t, null, null);
    } else
      o = I;
    return this.view.setOverlayImageLayer(o, E);
  }, B.prototype.getOverlayImageLayer = function(I = "overlay") {
    return this.view.getImageLayer(I);
  }, B.prototype.increaseZoom = function() {
    this.view.increaseZoom(0.01);
  }, B.prototype.decreaseZoom = function() {
    this.view.decreaseZoom(0.01);
  }, B.prototype.setRotation = function(I) {
    this.view.setRotation(I);
  }, B.prototype.setActiveHiPSLayer = function(I) {
    this.view.setActiveHiPSLayer(I);
  }, B.prototype.getActiveHiPSLayer = function() {
    return this.view.selectedLayer;
  }, B.prototype.getImageOverlays = function() {
    return this.view.overlayLayers;
  }, B.prototype.getOverlays = function() {
    return this.view.allOverlayLayers;
  }, B.prototype.getImageOverlays = function() {
    return this.view.overlayLayers;
  }, B.prototype.isHpxGridDisplayed = function() {
    return this.view.displayHpxGrid;
  }, B.prototype.isReticleDisplayed = function() {
    return this.view.displayReticle;
  }, B.prototype.createProgressiveCatalog = function(I, E, o, t) {
    return new NQ(I, E, o, t);
  }, B.prototype.createOverlay = function(I) {
    return new gI(I);
  }, B.AVAILABLE_CALLBACKS = [
    "select",
    "objectClicked",
    "objectHovered",
    "objectHoveredStop",
    "footprintClicked",
    "footprintHovered",
    "positionChanged",
    "zoomChanged",
    "click",
    "rightClickMove",
    "mouseMove",
    "fullScreenToggled",
    "cooFrameChanged"
  ], B.prototype.on = function(I, E) {
    B.AVAILABLE_CALLBACKS.indexOf(I) < 0 || (this.callbacksByEventName[I] = E, I === "positionChanged" && CA.AL_USE_WASM.dispatchedTo(document.body, { callback: (o) => {
      let t = Z.throttle(
        E,
        eI.CALLBACKS_THROTTLE_TIME_MS
      );
      o.setCallbackPositionChanged(t);
    } }));
  }, B.prototype.addListener = function(I, E) {
    new CA(I).listenedBy(this.aladinDiv, E);
  }, B.prototype.select = function() {
    this.fire("selectstart");
  }, B.prototype.fire = function(I, E) {
    if (I === "selectstart")
      this.view.setMode(eI.SELECT);
    else if (I === "selectend") {
      this.view.setMode(eI.PAN);
      var o = this.callbacksByEventName.select;
      typeof o == "function" && o(E);
    }
  }, B.prototype.hideBoxes = function() {
    if (this.boxes)
      for (var I = 0; I < this.boxes.length; I++)
        this.boxes[I].hide();
  }, B.prototype.updateCM = function() {
  }, B.prototype.showLayerBox = function() {
    this.stack.show();
  }, B.prototype.showCooGridBox = function() {
    this.coogrid.show();
  }, B.prototype.showCooGrid = function() {
    this.view.setGridConfig({ enabled: !0 }), this.cooGridEnabled = !0;
  }, B.prototype.hideCooGrid = function() {
    this.view.setGridConfig({ enabled: !1 }), this.cooGridEnabled = !1;
  }, B.prototype.layerByName = function(I) {
    for (var E = this.view.allOverlayLayers, o = 0; o < E.length; o++)
      if (I == E[o].name)
        return E[o];
    return null;
  }, B.prototype.exportAsPNG = function(I = !1) {
    (async () => {
      const E = await this.getViewDataURL();
      if (I)
        Z.download(E, "screenshot");
      else {
        var o = window.open();
        o.document.write('<img src="' + E + '" width="' + this.view.width + 'px">'), o.document.title = "Aladin Lite snapshot";
      }
    })();
  }, B.prototype.getViewDataURL = async function(E) {
    var E = E || {};
    if (typeof E != "object") {
      var o = E;
      E = { format: o };
    }
    return await this.view.getCanvasDataURL(E.format, E.width, E.height);
  }, B.prototype.getViewWCS = function(I) {
    var E = this.getRaDec(), o = this.getFov();
    return {
      NAXIS: 2,
      NAXIS1: this.view.width,
      NAXIS2: this.view.height,
      RADECSYS: "ICRS",
      CRPIX1: this.view.width / 2,
      CRPIX2: this.view.height / 2,
      CRVAL1: E[0],
      CRVAL2: E[1],
      CTYPE1: "RA---SIN",
      CTYPE2: "DEC--SIN",
      CD1_1: o[0] / this.view.width,
      CD1_2: 0,
      CD2_1: 0,
      CD2_2: o[1] / this.view.height
    };
  }, B.prototype.setFovRange = B.prototype.setFOVRange = function(I, E) {
    if (I > E) {
      var o = I;
      I = E, E = o;
    }
    this.view.minFOV = I, this.view.maxFOV = E;
  }, B.prototype.pix2world = function(I, E) {
    if (this.view)
      try {
        const [o, t] = this.wasm.screenToWorld(I, E);
        return o < 0 ? [o + 360, t] : [o, t];
      } catch {
        return;
      }
  }, B.prototype.world2pix = function(I, E) {
    if (this.view)
      try {
        return this.wasm.worldToScreen(I, E);
      } catch {
        return;
      }
  }, B.prototype.getFovCorners = function(I) {
    (!I || I < 1) && (I = 1);
    for (var E = [], o, t, s, M, N = 0; N < 4; N++) {
      o = N == 0 || N == 3 ? 0 : this.view.width - 1, t = N < 2 ? 0 : this.view.height - 1, s = N < 2 ? this.view.width - 1 : 0, M = N == 1 || N == 2 ? this.view.height - 1 : 0;
      for (var k = 0; k < I; k++) {
        let U = this.wasm.screenToWorld(o + k / I * (s - o), t + k / I * (M - t));
        E.push(U);
      }
    }
    return E;
  }, B.prototype.getFov = function() {
    var I = this.view.fov, E = this.getSize(), o = E[1] / E[0] * I;
    return I = Math.min(I, 180), o = Math.min(o, 180), [I, o];
  }, B.prototype.getSize = function() {
    return [this.view.width, this.view.height];
  }, B.prototype.getParentDiv = function() {
    return O(this.aladinDiv);
  }, B;
}();
VA.prototype.box = function(B) {
  var A = new Box(B);
  return A.$parentDiv.appendTo(this.aladinDiv), A;
};
VA.prototype.showPopup = function(B, A, g, C, I) {
  this.view.catalogForPopup.removeAll(), this.view.overlayForPopup.removeAll();
  let E;
  I !== void 0 ? (this.view.overlayForPopup.add(zA.circle(B, A, I, { fillColor: "rgba(255, 0, 0, 0.2)" })), E = zA.marker(B, A, { popupTitle: g, popupDesc: C, useMarkerDefaultIcon: !0 })) : E = zA.marker(B, A, { popupTitle: g, popupDesc: C, useMarkerDefaultIcon: !1 }), this.view.catalogForPopup.addSources(E), this.view.overlayForPopup.show(), this.view.catalogForPopup.show(), this.view.popup.setTitle(g), this.view.popup.setText(C), this.view.popup.setSource(E), this.view.popup.show();
};
VA.prototype.hidePopup = function() {
  this.view.popup.hide();
};
VA.prototype.getShareURL = function() {
  var B = this.getRaDec(), A = new $A();
  return A.prec = 7, A.lon = B[0], A.lat = B[1], VA.URL_PREVIEWER + "?target=" + encodeURIComponent(A.format("s")) + "&fov=" + this.getFov()[0].toFixed(2) + "&survey=" + encodeURIComponent(this.getBaseImageLayer().id || this.getBaseImageLayer().rootUrl);
};
VA.prototype.getEmbedCode = function() {
  var B = this.getRaDec(), A = new $A();
  A.prec = 7, A.lon = B[0], A.lat = B[1], this.getBaseImageLayer().id;
  var g = this.getFov()[0];
  let C = "";
  const I = `
`;
  return C += '<div id="aladin-lite-div" style="width:400px;height:400px;"></div>' + I, C += '<script src="https://aladin.cds.unistra.fr/AladinLite/api/v3/latest/aladin.js" charset="utf-8"><\/script>' + I, C += "<script>" + I, C += "let aladin;" + I + "A.init.then(() => {" + I + "   aladin = A.aladin('#aladin-lite-div', {survey: 'P/DSS2/color', fov: " + g.toFixed(2) + ', target: "' + A.format("s") + '"});' + I + "});" + I, C += "<\/script>", C;
};
VA.prototype.displayFITS = function(B, A, g, C, I = "base") {
  g = g || ((o, t, s, M) => {
    this.gotoRaDec(o, t), this.setFoV(s);
  });
  const E = this.createImageFITS(B, B, A, g, C);
  return this.setOverlayImageLayer(E, I);
};
VA.prototype.displayJPG = VA.prototype.displayPNG = function(B, A, g, C) {
  A = A || {}, A.color = !0, A.label = "JPG/PNG image", A.outputFormat = "png", A = A || {};
  var I = { url: B };
  A.color && (I.color = !0), A.outputFormat && (I.format = A.outputFormat), A.order && (I.order = A.order), A.nocache && (I.nocache = A.nocache);
  let E = this;
  const o = (s, M = {}, N = "GET") => {
    let k = {
      method: N
    };
    return N === "GET" ? s += "?" + new URLSearchParams(M).toString() : k.body = JSON.stringify(M), fetch(s, k).then((U) => U.json());
  };
  ((s, M) => o(s, M, "GET"))("https://alasky.unistra.fr/cgi/fits2HiPS", I).then(async (s) => {
    if (s.status != "success") {
      console.error("An error occured: " + s.message), C && C(s.message);
      return;
    }
    var M = A.label || "FITS image", N = s.data.meta;
    const k = E.createImageSurvey(s.data.url, M, s.data.url);
    E.setOverlayImageLayer(k, "overlay"), A && A.transparency;
    var U = !0;
    g && (U = g(N.ra, N.dec, N.fov)), U === !0 && (E.wasm.setCenter(N.ra, N.dec), E.setFoV(N.fov));
  });
};
VA.prototype.setReduceDeformations = function(B) {
  this.reduceDeformations = B, this.view.requestRedraw();
};
let RQ = function() {
  let B = function(A) {
    this.uuid = Z.uuidv4(), this.type = "moc", A = A || {}, this.name = A.name || "MOC", this.color = A.color || _A.getNextColor(), this.color = _A.standardizeColor(this.color), this.fillColor = A.fillColor || this.color, this.fillColor = _A.standardizeColor(this.fillColor), A && A.perimeter ? this.perimeter = !0 : this.perimeter = !1, A && A.fill ? this.fill = !0 : this.fill = !1, A && A.edge ? this.edge = !0 : this.edge = !1, !this.fill && !this.perimeter && A && !A.edge && (this.edge = !0), this.opacity = A.opacity || 1, this.opacity = Math.max(0, Math.min(1, this.opacity)), this.lineWidth = A.lineWidth || 1, this.isShowing = !0, this.ready = !1, this.skyFrac = void 0;
  };
  return B.prototype.skyFraction = function() {
    return this.skyFrac;
  }, B.prototype.dataFromJSON = function(A) {
    this.dataJSON = A;
  }, B.prototype.dataFromFITSURL = function(A, g) {
    this.dataURL = A, this.promiseFetchData = fetch(this.dataURL).then((C) => C.arrayBuffer()), this.successCallback = g;
  }, B.prototype.setView = function(A) {
    let g = this;
    this.view = A, this.mocParams = new VA.wasmLibs.core.MOC(this.uuid, this.opacity, this.lineWidth, this.perimeter, this.fill, this.edge, this.isShowing, this.color, this.fillColor), this.dataURL ? this.promiseFetchData.then((C) => {
      g.view.wasm.addFITSMoc(g.mocParams, new Uint8Array(C)), g.ready = !0, g.successCallback && g.successCallback(g), g.skyFrac = g.view.wasm.mocSkyFraction(this.mocParams), g.view.mocs.push(g), g.view.allOverlayLayers.push(g), CA.GRAPHIC_OVERLAY_LAYER_ADDED.dispatchedTo(g.view.aladinDiv, { layer: g }), g.view.requestRedraw();
    }) : this.dataFromJSON && (g.view.wasm.addJSONMoc(g.mocParams, g.dataJSON), g.ready = !0, g.skyFrac = g.view.wasm.mocSkyFraction(g.mocParams), g.view.mocs.push(g), g.view.allOverlayLayers.push(g), CA.GRAPHIC_OVERLAY_LAYER_ADDED.dispatchedTo(g.view.aladinDiv, { layer: g }), g.view.requestRedraw());
  }, B.prototype.reportChange = function() {
    this.view && (this.mocParams = new VA.wasmLibs.core.MOC(this.uuid, this.opacity, this.lineWidth, this.perimeter, this.fill, this.edge, this.isShowing, this.color, this.fillColor), this.view.wasm.setMocParams(this.mocParams), this.view.requestRedraw());
  }, B.prototype.delete = function() {
    this.view && (this.view.wasm.removeMoc(this.mocParams), this.view.requestRedraw());
  }, B.prototype.show = function() {
    this.isShowing || (this.isShowing = !0, this.reportChange());
  }, B.prototype.hide = function() {
    this.isShowing && (this.isShowing = !1, this.reportChange());
  }, B.prototype.contains = function(A, g) {
    if (!this.ready)
      throw this.name + " is not yet ready, either because it has not been downloaded yet or because it has not been added to the aladin instance.";
    return this.view.wasm.mocContains(this.mocParams, A, g);
  }, B;
}(), zg = function() {
  return {
    buildSimbadCSURL: function(A, g) {
      if (A && typeof A == "object" && "ra" in A && "dec" in A) {
        var C = new $A(A.ra, A.dec, 7);
        A = C.format("s");
      }
      return "https://alasky.unistra.fr/cgi/simbad-flat/simbad-cs.py?target=" + encodeURIComponent(A) + "&SR=" + g + "&format=votable&SRUNIT=deg&SORTBY=nbref";
    },
    buildNEDPositionCSURL: function(A, g, C) {
      return "https://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?search_type=Near+Position+Search&of=xml_main&RA=" + A + "&DEC=" + g + "&SR=" + C;
    },
    buildNEDObjectCSURL: function(A, g) {
      return "https://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?search_type=Near+Name+Search&radius=" + 60 * g + "&of=xml_main&objname=" + A;
    },
    buildVizieRCSURL: function(A, g, C, I) {
      if (g && typeof g == "object" && "ra" in g && "dec" in g) {
        var E = new $A(g.ra, g.dec, 7);
        g = E.format("s");
      }
      var o = 1e5;
      I && I.hasOwnProperty("limit") && Z.isNumber(I.limit) && (o = parseInt(I.limit));
      let t = "https://vizier.unistra.fr/viz-bin/votable?-source=" + A + "&-c=" + encodeURIComponent(g) + "&-out.max=" + o + "&-c.rd=" + C;
      return t = t + "&-out.add=s_region", t = t + "&-out.add=s_fov", t;
    },
    buildSkyBotCSURL: function(A, g, C, I, E) {
      var o = "http://vo.imcce.fr/webservices/skybot/skybotconesearch_query.php?-from=AladinLite";
      if (o += "&RA=" + encodeURIComponent(A), o += "&DEC=" + encodeURIComponent(g), o += "&SR=" + encodeURIComponent(C), o += "&EPOCH=" + encodeURIComponent(I), E)
        for (var t in E)
          E.hasOwnProperty(t) && (o += "&" + t + "=" + encodeURIComponent(E[t]));
      return o;
    }
  };
}(), h;
const KQ = typeof TextDecoder < "u" ? new TextDecoder("utf-8", { ignoreBOM: !0, fatal: !0 }) : { decode: () => {
  throw Error("TextDecoder not available");
} };
typeof TextDecoder < "u" && KQ.decode();
let cI = null;
function $g() {
  return (cI === null || cI.byteLength === 0) && (cI = new Uint8Array(h.memory.buffer)), cI;
}
function dA(B, A) {
  return B = B >>> 0, KQ.decode($g().subarray(B, B + A));
}
const Mg = new Array(128).fill(void 0);
Mg.push(void 0, null, !0, !1);
let rI = Mg.length;
function x(B) {
  rI === Mg.length && Mg.push(Mg.length + 1);
  const A = rI;
  return rI = Mg[A], Mg[A] = B, A;
}
function r(B) {
  return Mg[B];
}
function uE(B) {
  B < 132 || (Mg[B] = rI, rI = B);
}
function AA(B) {
  const A = r(B);
  return uE(B), A;
}
function FA(B) {
  return B == null;
}
let GI = null;
function AI() {
  return (GI === null || GI.byteLength === 0) && (GI = new Float64Array(h.memory.buffer)), GI;
}
let hI = null;
function f() {
  return (hI === null || hI.byteLength === 0) && (hI = new Int32Array(h.memory.buffer)), hI;
}
let MA = 0;
const xI = typeof TextEncoder < "u" ? new TextEncoder("utf-8") : { encode: () => {
  throw Error("TextEncoder not available");
} }, xE = typeof xI.encodeInto == "function" ? function(B, A) {
  return xI.encodeInto(B, A);
} : function(B, A) {
  const g = xI.encode(B);
  return A.set(g), {
    read: B.length,
    written: g.length
  };
};
function JA(B, A, g) {
  if (g === void 0) {
    const t = xI.encode(B), s = A(t.length, 1) >>> 0;
    return $g().subarray(s, s + t.length).set(t), MA = t.length, s;
  }
  let C = B.length, I = A(C, 1) >>> 0;
  const E = $g();
  let o = 0;
  for (; o < C; o++) {
    const t = B.charCodeAt(o);
    if (t > 127)
      break;
    E[I + o] = t;
  }
  if (o !== C) {
    o !== 0 && (B = B.slice(o)), I = g(I, C, C = o + B.length * 3, 1) >>> 0;
    const t = $g().subarray(I + o, I + C), s = xE(B, t);
    o += s.written;
  }
  return MA = o, I;
}
function sB(B) {
  const A = typeof B;
  if (A == "number" || A == "boolean" || B == null)
    return `${B}`;
  if (A == "string")
    return `"${B}"`;
  if (A == "symbol") {
    const I = B.description;
    return I == null ? "Symbol" : `Symbol(${I})`;
  }
  if (A == "function") {
    const I = B.name;
    return typeof I == "string" && I.length > 0 ? `Function(${I})` : "Function";
  }
  if (Array.isArray(B)) {
    const I = B.length;
    let E = "[";
    I > 0 && (E += sB(B[0]));
    for (let o = 1; o < I; o++)
      E += ", " + sB(B[o]);
    return E += "]", E;
  }
  const g = /\[object ([^\]]+)\]/.exec(toString.call(B));
  let C;
  if (g.length > 1)
    C = g[1];
  else
    return toString.call(B);
  if (C == "Object")
    try {
      return "Object(" + JSON.stringify(B) + ")";
    } catch {
      return "Object";
    }
  return B instanceof Error ? `${B.name}: ${B.message}
${B.stack}` : C;
}
function oQ(B, A, g, C) {
  const I = { a: B, b: A, cnt: 1, dtor: g }, E = (...o) => {
    I.cnt++;
    const t = I.a;
    I.a = 0;
    try {
      return C(t, I.b, ...o);
    } finally {
      --I.cnt === 0 ? h.__wbindgen_export_2.get(I.dtor)(t, I.b) : I.a = t;
    }
  };
  return E.original = I, E;
}
function mE(B, A) {
  h.wasm_bindgen__convert__closures__invoke0_mut__hfb7d646e97c7b2fa(B, A);
}
function OE(B, A, g) {
  h.wasm_bindgen__convert__closures__invoke1_mut__h1c5403ccd99fb092(B, A, x(g));
}
function vE(B, A, g, C) {
  const I = { a: B, b: A, cnt: 1, dtor: g }, E = (...o) => {
    I.cnt++;
    try {
      return C(I.a, I.b, ...o);
    } finally {
      --I.cnt === 0 && (h.__wbindgen_export_2.get(I.dtor)(I.a, I.b), I.a = 0);
    }
  };
  return E.original = I, E;
}
function bE(B, A) {
  h._dyn_core__ops__function__Fn_____Output___R_as_wasm_bindgen__closure__WasmClosure___describe__invoke__hbd7155e487187523(B, A);
}
function Pg(B, A) {
  return B = B >>> 0, AI().subarray(B / 8, B / 8 + A);
}
function DQ(B, A) {
  const g = A(B.length * 8, 8) >>> 0;
  return AI().set(B, g / 8), MA = B.length, g;
}
let NI = 128;
function aQ(B) {
  if (NI == 1)
    throw new Error("out of js stack");
  return Mg[--NI] = B, NI;
}
let MI = null;
function YQ() {
  return (MI === null || MI.byteLength === 0) && (MI = new Uint32Array(h.memory.buffer)), MI;
}
function ZE(B, A) {
  B = B >>> 0;
  const C = YQ().subarray(B / 4, B / 4 + A), I = [];
  for (let E = 0; E < C.length; E++)
    I.push(AA(C[E]));
  return I;
}
function TE(B, A) {
  const g = A(B.length * 4, 4) >>> 0, C = YQ();
  for (let I = 0; I < B.length; I++)
    C[g / 4 + I] = x(B[I]);
  return MA = B.length, g;
}
function kg(B, A) {
  if (!(B instanceof A))
    throw new Error(`expected instance of ${A.name}`);
  return B.ptr;
}
function lQ(B, A) {
  const g = A(B.length * 1, 1) >>> 0;
  return $g().set(B, g / 1), MA = B.length, g;
}
let yI = null;
function kI() {
  return (yI === null || yI.byteLength === 0) && (yI = new Float32Array(h.memory.buffer)), yI;
}
function kA(B, A) {
  try {
    return B.apply(this, A);
  } catch (g) {
    h.__wbindgen_exn_store(x(g));
  }
}
function WE(B, A, g, C) {
  h.wasm_bindgen__convert__closures__invoke2_mut__h76c8ca53f22afe4e(B, A, x(g), x(C));
}
function VE(B, A) {
  return B = B >>> 0, $g().subarray(B / 1, B / 1 + A);
}
function jE(B, A) {
  return B = B >>> 0, kI().subarray(B / 4, B / 4 + A);
}
const _E = Object.freeze({ Fits: 0, 0: "Fits", Jpeg: 1, 1: "Jpeg", Png: 2, 2: "Png", Webp: 3, 3: "Webp" }), PE = Object.freeze({ Linear: 0, 0: "Linear", Sqrt: 1, 1: "Sqrt", Log: 2, 2: "Log", Asinh: 3, 3: "Asinh", Pow2: 4, 4: "Pow2" }), XE = Object.freeze({ DMM: 0, 0: "DMM", DD: 1, 1: "DD", DMS: 2, 2: "DMS", HMS: 3, 3: "HMS" }), zE = Object.freeze({ ICRS: 0, 0: "ICRS", GAL: 1, 1: "GAL" }), $E = Object.freeze({ Zero: 0, 0: "Zero", One: 1, 1: "One", SrcColor: 2, 2: "SrcColor", OneMinusSrcColor: 3, 3: "OneMinusSrcColor", DstColor: 4, 4: "DstColor", OneMinusDstColor: 5, 5: "OneMinusDstColor", SrcAlpha: 6, 6: "SrcAlpha", OneMinusSrcAlpha: 7, 7: "OneMinusSrcAlpha", DstAlpha: 8, 8: "DstAlpha", OneMinusDstAlpha: 9, 9: "OneMinusDstAlpha", ConstantColor: 10, 10: "ConstantColor", OneMinusConstantColor: 11, 11: "OneMinusConstantColor", ConstantAlpha: 12, 12: "ConstantAlpha", OneMinusConstantAlpha: 13, 13: "OneMinusConstantAlpha" }), Ai = Object.freeze({ FuncAdd: 0, 0: "FuncAdd", FuncSubstract: 1, 1: "FuncSubstract", FuncReverseSubstract: 2, 2: "FuncReverseSubstract" });
class nI {
  static __wrap(A) {
    A = A >>> 0;
    const g = Object.create(nI.prototype);
    return g.__wbg_ptr = A, g;
  }
  __destroy_into_raw() {
    const A = this.__wbg_ptr;
    return this.__wbg_ptr = 0, A;
  }
  free() {
    const A = this.__destroy_into_raw();
    h.__wbg_blendcfg_free(A);
  }
  /**
  * @returns {number}
  */
  get src_color_factor() {
    return h.__wbg_get_blendcfg_src_color_factor(this.__wbg_ptr) >>> 0;
  }
  /**
  * @param {number} arg0
  */
  set src_color_factor(A) {
    h.__wbg_set_blendcfg_src_color_factor(this.__wbg_ptr, A);
  }
  /**
  * @returns {number}
  */
  get dst_color_factor() {
    return h.__wbg_get_blendcfg_dst_color_factor(this.__wbg_ptr) >>> 0;
  }
  /**
  * @param {number} arg0
  */
  set dst_color_factor(A) {
    h.__wbg_set_blendcfg_dst_color_factor(this.__wbg_ptr, A);
  }
  /**
  * @returns {number}
  */
  get func() {
    return h.__wbg_get_blendcfg_func(this.__wbg_ptr) >>> 0;
  }
  /**
  * @param {number} arg0
  */
  set func(A) {
    h.__wbg_set_blendcfg_func(this.__wbg_ptr, A);
  }
}
class gi {
  __destroy_into_raw() {
    const A = this.__wbg_ptr;
    return this.__wbg_ptr = 0, A;
  }
  free() {
    const A = this.__destroy_into_raw();
    h.__wbg_centeredfov_free(A);
  }
  /**
  * Position of the field of view
  * @returns {number}
  */
  get ra() {
    return h.__wbg_get_centeredfov_ra(this.__wbg_ptr);
  }
  /**
  * Position of the field of view
  * @param {number} arg0
  */
  set ra(A) {
    h.__wbg_set_centeredfov_ra(this.__wbg_ptr, A);
  }
  /**
  * @returns {number}
  */
  get dec() {
    return h.__wbg_get_centeredfov_dec(this.__wbg_ptr);
  }
  /**
  * @param {number} arg0
  */
  set dec(A) {
    h.__wbg_set_centeredfov_dec(this.__wbg_ptr, A);
  }
  /**
  * Aperture
  * @returns {number}
  */
  get fov() {
    return h.__wbg_get_centeredfov_fov(this.__wbg_ptr);
  }
  /**
  * Aperture
  * @param {number} arg0
  */
  set fov(A) {
    h.__wbg_set_centeredfov_fov(this.__wbg_ptr, A);
  }
}
class FI {
  static __wrap(A) {
    A = A >>> 0;
    const g = Object.create(FI.prototype);
    return g.__wbg_ptr = A, g;
  }
  __destroy_into_raw() {
    const A = this.__wbg_ptr;
    return this.__wbg_ptr = 0, A;
  }
  free() {
    const A = this.__destroy_into_raw();
    h.__wbg_colorrgb_free(A);
  }
  /**
  * @returns {number}
  */
  get r() {
    return h.__wbg_get_colorrgb_r(this.__wbg_ptr);
  }
  /**
  * @param {number} arg0
  */
  set r(A) {
    h.__wbg_set_colorrgb_r(this.__wbg_ptr, A);
  }
  /**
  * @returns {number}
  */
  get g() {
    return h.__wbg_get_colorrgb_g(this.__wbg_ptr);
  }
  /**
  * @param {number} arg0
  */
  set g(A) {
    h.__wbg_set_colorrgb_g(this.__wbg_ptr, A);
  }
  /**
  * @returns {number}
  */
  get b() {
    return h.__wbg_get_colorrgb_b(this.__wbg_ptr);
  }
  /**
  * @param {number} arg0
  */
  set b(A) {
    h.__wbg_set_colorrgb_b(this.__wbg_ptr, A);
  }
}
class vg {
  static __wrap(A) {
    A = A >>> 0;
    const g = Object.create(vg.prototype);
    return g.__wbg_ptr = A, g;
  }
  __destroy_into_raw() {
    const A = this.__wbg_ptr;
    return this.__wbg_ptr = 0, A;
  }
  free() {
    const A = this.__destroy_into_raw();
    h.__wbg_colorrgba_free(A);
  }
  /**
  * @returns {number}
  */
  get r() {
    return h.__wbg_get_colorrgb_r(this.__wbg_ptr);
  }
  /**
  * @param {number} arg0
  */
  set r(A) {
    h.__wbg_set_colorrgb_r(this.__wbg_ptr, A);
  }
  /**
  * @returns {number}
  */
  get g() {
    return h.__wbg_get_colorrgb_g(this.__wbg_ptr);
  }
  /**
  * @param {number} arg0
  */
  set g(A) {
    h.__wbg_set_colorrgb_g(this.__wbg_ptr, A);
  }
  /**
  * @returns {number}
  */
  get b() {
    return h.__wbg_get_colorrgb_b(this.__wbg_ptr);
  }
  /**
  * @param {number} arg0
  */
  set b(A) {
    h.__wbg_set_colorrgb_b(this.__wbg_ptr, A);
  }
  /**
  * @returns {number}
  */
  get a() {
    return h.__wbg_get_colorrgba_a(this.__wbg_ptr);
  }
  /**
  * @param {number} arg0
  */
  set a(A) {
    h.__wbg_set_colorrgba_a(this.__wbg_ptr, A);
  }
}
class Ii {
  __destroy_into_raw() {
    const A = this.__wbg_ptr;
    return this.__wbg_ptr = 0, A;
  }
  free() {
    const A = this.__destroy_into_raw();
    h.__wbg_gridcfg_free(A);
  }
  /**
  * @returns {ColorRGB | undefined}
  */
  get color() {
    const A = h.__wbg_get_gridcfg_color(this.__wbg_ptr);
    return A === 0 ? void 0 : FI.__wrap(A);
  }
  /**
  * @param {ColorRGB | undefined} arg0
  */
  set color(A) {
    let g = 0;
    FA(A) || (kg(A, FI), g = A.__destroy_into_raw()), h.__wbg_set_gridcfg_color(this.__wbg_ptr, g);
  }
  /**
  * @returns {number | undefined}
  */
  get thickness() {
    try {
      const C = h.__wbindgen_add_to_stack_pointer(-16);
      h.__wbg_get_gridcfg_thickness(C, this.__wbg_ptr);
      var A = f()[C / 4 + 0], g = kI()[C / 4 + 1];
      return A === 0 ? void 0 : g;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {number | undefined} arg0
  */
  set thickness(A) {
    h.__wbg_set_gridcfg_thickness(this.__wbg_ptr, !FA(A), FA(A) ? 0 : A);
  }
  /**
  * @returns {number | undefined}
  */
  get opacity() {
    try {
      const C = h.__wbindgen_add_to_stack_pointer(-16);
      h.__wbg_get_gridcfg_opacity(C, this.__wbg_ptr);
      var A = f()[C / 4 + 0], g = kI()[C / 4 + 1];
      return A === 0 ? void 0 : g;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {number | undefined} arg0
  */
  set opacity(A) {
    h.__wbg_set_gridcfg_opacity(this.__wbg_ptr, !FA(A), FA(A) ? 0 : A);
  }
  /**
  * @returns {boolean | undefined}
  */
  get show_labels() {
    const A = h.__wbg_get_gridcfg_show_labels(this.__wbg_ptr);
    return A === 16777215 ? void 0 : A !== 0;
  }
  /**
  * @param {boolean | undefined} arg0
  */
  set show_labels(A) {
    h.__wbg_set_gridcfg_show_labels(this.__wbg_ptr, FA(A) ? 16777215 : A ? 1 : 0);
  }
  /**
  * @returns {number | undefined}
  */
  get label_size() {
    try {
      const C = h.__wbindgen_add_to_stack_pointer(-16);
      h.__wbg_get_gridcfg_label_size(C, this.__wbg_ptr);
      var A = f()[C / 4 + 0], g = kI()[C / 4 + 1];
      return A === 0 ? void 0 : g;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {number | undefined} arg0
  */
  set label_size(A) {
    h.__wbg_set_gridcfg_label_size(this.__wbg_ptr, !FA(A), FA(A) ? 0 : A);
  }
  /**
  * @returns {boolean | undefined}
  */
  get enabled() {
    const A = h.__wbg_get_gridcfg_enabled(this.__wbg_ptr);
    return A === 16777215 ? void 0 : A !== 0;
  }
  /**
  * @param {boolean | undefined} arg0
  */
  set enabled(A) {
    h.__wbg_set_gridcfg_enabled(this.__wbg_ptr, FA(A) ? 16777215 : A ? 1 : 0);
  }
  /**
  * @returns {number | undefined}
  */
  get fmt() {
    const A = h.__wbg_get_gridcfg_fmt(this.__wbg_ptr);
    return A === 4 ? void 0 : A;
  }
  /**
  * @param {number | undefined} arg0
  */
  set fmt(A) {
    h.__wbg_set_gridcfg_fmt(this.__wbg_ptr, FA(A) ? 4 : A);
  }
}
class vI {
  static __wrap(A) {
    A = A >>> 0;
    const g = Object.create(vI.prototype);
    return g.__wbg_ptr = A, g;
  }
  __destroy_into_raw() {
    const A = this.__wbg_ptr;
    return this.__wbg_ptr = 0, A;
  }
  free() {
    const A = this.__destroy_into_raw();
    h.__wbg_imagemetadata_free(A);
  }
  /**
  * @returns {BlendCfg}
  */
  get blend_cfg() {
    const A = h.__wbg_get_imagemetadata_blend_cfg(this.__wbg_ptr);
    return nI.__wrap(A);
  }
  /**
  * @param {BlendCfg} arg0
  */
  set blend_cfg(A) {
    kg(A, nI);
    var g = A.__destroy_into_raw();
    h.__wbg_set_imagemetadata_blend_cfg(this.__wbg_ptr, g);
  }
  /**
  * @returns {number}
  */
  get opacity() {
    return h.__wbg_get_imagemetadata_opacity(this.__wbg_ptr);
  }
  /**
  * @param {number} arg0
  */
  set opacity(A) {
    h.__wbg_set_imagemetadata_opacity(this.__wbg_ptr, A);
  }
  /**
  * @returns {boolean}
  */
  get longitude_reversed() {
    return h.__wbg_get_imagemetadata_longitude_reversed(this.__wbg_ptr) !== 0;
  }
  /**
  * @param {boolean} arg0
  */
  set longitude_reversed(A) {
    h.__wbg_set_imagemetadata_longitude_reversed(this.__wbg_ptr, A);
  }
  /**
  * the current format chosen
  * @returns {number}
  */
  get img_format() {
    return h.__wbg_get_imagemetadata_img_format(this.__wbg_ptr) >>> 0;
  }
  /**
  * the current format chosen
  * @param {number} arg0
  */
  set img_format(A) {
    h.__wbg_set_imagemetadata_img_format(this.__wbg_ptr, A);
  }
  /**
  * @param {any} color
  */
  set color(A) {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16);
      h.imagemetadata_set_color(I, this.__wbg_ptr, x(A));
      var g = f()[I / 4 + 0], C = f()[I / 4 + 1];
      if (C)
        throw AA(g);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @returns {any}
  */
  get color() {
    const A = h.imagemetadata_color(this.__wbg_ptr);
    return AA(A);
  }
}
class Bi {
  __destroy_into_raw() {
    const A = this.__wbg_ptr;
    return this.__wbg_ptr = 0, A;
  }
  free() {
    const A = this.__destroy_into_raw();
    h.__wbg_intounderlyingbytesource_free(A);
  }
  /**
  * @returns {string}
  */
  get type() {
    let A, g;
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16);
      h.intounderlyingbytesource_type(E, this.__wbg_ptr);
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      return A = C, g = I, dA(C, I);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16), h.__wbindgen_free(A, g, 1);
    }
  }
  /**
  * @returns {number}
  */
  get autoAllocateChunkSize() {
    return h.intounderlyingbytesource_autoAllocateChunkSize(this.__wbg_ptr) >>> 0;
  }
  /**
  * @param {any} controller
  */
  start(A) {
    h.intounderlyingbytesource_start(this.__wbg_ptr, x(A));
  }
  /**
  * @param {any} controller
  * @returns {Promise<any>}
  */
  pull(A) {
    const g = h.intounderlyingbytesource_pull(this.__wbg_ptr, x(A));
    return AA(g);
  }
  /**
  */
  cancel() {
    const A = this.__destroy_into_raw();
    h.__wbg_intounderlyingbytesource_free(A);
  }
}
class Qi {
  __destroy_into_raw() {
    const A = this.__wbg_ptr;
    return this.__wbg_ptr = 0, A;
  }
  free() {
    const A = this.__destroy_into_raw();
    h.__wbg_intounderlyingsink_free(A);
  }
  /**
  * @param {any} chunk
  * @returns {Promise<any>}
  */
  write(A) {
    const g = h.intounderlyingsink_write(this.__wbg_ptr, x(A));
    return AA(g);
  }
  /**
  * @returns {Promise<any>}
  */
  close() {
    const A = this.__destroy_into_raw(), g = h.intounderlyingsink_close(A);
    return AA(g);
  }
  /**
  * @param {any} reason
  * @returns {Promise<any>}
  */
  abort(A) {
    const g = this.__destroy_into_raw(), C = h.intounderlyingsink_abort(g, x(A));
    return AA(C);
  }
}
class Ci {
  __destroy_into_raw() {
    const A = this.__wbg_ptr;
    return this.__wbg_ptr = 0, A;
  }
  free() {
    const A = this.__destroy_into_raw();
    h.__wbg_intounderlyingsource_free(A);
  }
  /**
  * @param {any} controller
  * @returns {Promise<any>}
  */
  pull(A) {
    const g = h.intounderlyingsource_pull(this.__wbg_ptr, x(A));
    return AA(g);
  }
  /**
  */
  cancel() {
    const A = this.__destroy_into_raw();
    h.intounderlyingsource_cancel(A);
  }
}
class ng {
  static __wrap(A) {
    A = A >>> 0;
    const g = Object.create(ng.prototype);
    return g.__wbg_ptr = A, g;
  }
  __destroy_into_raw() {
    const A = this.__wbg_ptr;
    return this.__wbg_ptr = 0, A;
  }
  free() {
    const A = this.__destroy_into_raw();
    h.__wbg_moc_free(A);
  }
  /**
  * @returns {number}
  */
  get line_width() {
    return h.__wbg_get_moc_line_width(this.__wbg_ptr);
  }
  /**
  * @param {number} arg0
  */
  set line_width(A) {
    h.__wbg_set_moc_line_width(this.__wbg_ptr, A);
  }
  /**
  * @returns {boolean}
  */
  get perimeter() {
    return h.__wbg_get_moc_perimeter(this.__wbg_ptr) !== 0;
  }
  /**
  * @param {boolean} arg0
  */
  set perimeter(A) {
    h.__wbg_set_moc_perimeter(this.__wbg_ptr, A);
  }
  /**
  * @returns {boolean}
  */
  get filled() {
    return h.__wbg_get_moc_filled(this.__wbg_ptr) !== 0;
  }
  /**
  * @param {boolean} arg0
  */
  set filled(A) {
    h.__wbg_set_moc_filled(this.__wbg_ptr, A);
  }
  /**
  * @returns {boolean}
  */
  get edges() {
    return h.__wbg_get_moc_edges(this.__wbg_ptr) !== 0;
  }
  /**
  * @param {boolean} arg0
  */
  set edges(A) {
    h.__wbg_set_moc_edges(this.__wbg_ptr, A);
  }
  /**
  * @returns {boolean}
  */
  get show() {
    return h.__wbg_get_moc_show(this.__wbg_ptr) !== 0;
  }
  /**
  * @param {boolean} arg0
  */
  set show(A) {
    h.__wbg_set_moc_show(this.__wbg_ptr, A);
  }
  /**
  * @returns {ColorRGBA}
  */
  get color() {
    const A = h.__wbg_get_moc_color(this.__wbg_ptr);
    return vg.__wrap(A);
  }
  /**
  * @param {ColorRGBA} arg0
  */
  set color(A) {
    kg(A, vg);
    var g = A.__destroy_into_raw();
    h.__wbg_set_moc_color(this.__wbg_ptr, g);
  }
  /**
  * @returns {ColorRGBA}
  */
  get fill_color() {
    const A = h.__wbg_get_moc_fill_color(this.__wbg_ptr);
    return vg.__wrap(A);
  }
  /**
  * @param {ColorRGBA} arg0
  */
  set fill_color(A) {
    kg(A, vg);
    var g = A.__destroy_into_raw();
    h.__wbg_set_moc_fill_color(this.__wbg_ptr, g);
  }
  /**
  * @param {string} uuid
  * @param {number} opacity
  * @param {number} line_width
  * @param {boolean} perimeter
  * @param {boolean} filled
  * @param {boolean} edges
  * @param {boolean} show
  * @param {string} hex_color
  * @param {string} fill_color
  */
  constructor(A, g, C, I, E, o, t, s, M) {
    const N = JA(A, h.__wbindgen_malloc, h.__wbindgen_realloc), k = MA, U = JA(s, h.__wbindgen_malloc, h.__wbindgen_realloc), Y = MA, K = JA(M, h.__wbindgen_malloc, h.__wbindgen_realloc), l = MA, d = h.moc_new(N, k, g, C, I, E, o, t, U, Y, K, l);
    return ng.__wrap(d);
  }
}
class Ei {
  __destroy_into_raw() {
    const A = this.__wbg_ptr;
    return this.__wbg_ptr = 0, A;
  }
  free() {
    const A = this.__destroy_into_raw();
    h.__wbg_pipeoptions_free(A);
  }
  /**
  * @returns {boolean}
  */
  get preventClose() {
    return h.pipeoptions_preventClose(this.__wbg_ptr) !== 0;
  }
  /**
  * @returns {boolean}
  */
  get preventCancel() {
    return h.pipeoptions_preventCancel(this.__wbg_ptr) !== 0;
  }
  /**
  * @returns {boolean}
  */
  get preventAbort() {
    return h.pipeoptions_preventAbort(this.__wbg_ptr) !== 0;
  }
  /**
  * @returns {AbortSignal | undefined}
  */
  get signal() {
    const A = h.pipeoptions_signal(this.__wbg_ptr);
    return AA(A);
  }
}
class ii {
  __destroy_into_raw() {
    const A = this.__wbg_ptr;
    return this.__wbg_ptr = 0, A;
  }
  free() {
    const A = this.__destroy_into_raw();
    h.__wbg_queuingstrategy_free(A);
  }
  /**
  * @returns {number}
  */
  get highWaterMark() {
    return h.queuingstrategy_highWaterMark(this.__wbg_ptr);
  }
}
class bI {
  static __wrap(A) {
    A = A >>> 0;
    const g = Object.create(bI.prototype);
    return g.__wbg_ptr = A, g;
  }
  __destroy_into_raw() {
    const A = this.__wbg_ptr;
    return this.__wbg_ptr = 0, A;
  }
  free() {
    const A = this.__destroy_into_raw();
    h.__wbg_readablestreamgetreaderoptions_free(A);
  }
  /**
  * @returns {any}
  */
  get mode() {
    const A = h.readablestreamgetreaderoptions_mode(this.__wbg_ptr);
    return AA(A);
  }
}
class mI {
  static __wrap(A) {
    A = A >>> 0;
    const g = Object.create(mI.prototype);
    return g.__wbg_ptr = A, g;
  }
  __destroy_into_raw() {
    const A = this.__wbg_ptr;
    return this.__wbg_ptr = 0, A;
  }
  free() {
    const A = this.__destroy_into_raw();
    h.__wbg_webclient_free(A);
  }
  /**
  * Create the Aladin Lite webgl backend
  *
  * # Arguments
  *
  * * `aladin_div_name` - The name of the div where aladin is created
  * * `shaders` - The list of shader objects containing the GLSL code source
  * * `resources` - Image resource files
  * @param {string} aladin_div_name
  * @param {any} shaders
  * @param {any} resources
  */
  constructor(A, g, C) {
    try {
      const t = h.__wbindgen_add_to_stack_pointer(-16), s = JA(A, h.__wbindgen_malloc, h.__wbindgen_realloc), M = MA;
      h.webclient_new(t, s, M, x(g), x(C));
      var I = f()[t / 4 + 0], E = f()[t / 4 + 1], o = f()[t / 4 + 2];
      if (o)
        throw AA(E);
      return mI.__wrap(I);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {Function} callback
  */
  setCallbackPositionChanged(A) {
    h.webclient_setCallbackPositionChanged(this.__wbg_ptr, x(A));
  }
  /**
  * Update the view
  *
  * # Arguments
  *
  * * `dt` - The time elapsed from the last frame update
  * * `force` - This parameter ensures to force the update of some elements
  *   even if the camera has not moved
  *
  * # Return
  * Whether the view is moving or not
  * @param {number} dt
  * @returns {boolean}
  */
  update(A) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_update(E, this.__wbg_ptr, A);
      var g = f()[E / 4 + 0], C = f()[E / 4 + 1], I = f()[E / 4 + 2];
      if (I)
        throw AA(C);
      return g !== 0;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Resize the window
  *
  * # Arguments
  *
  * * `width` - The width in pixels of the view
  * * `height` - The height in pixels of the view
  * @param {number} width
  * @param {number} height
  */
  resize(A, g) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_resize(E, this.__wbg_ptr, A, g);
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      if (I)
        throw AA(C);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Set the type of projections
  *
  * # Arguments
  *
  * * `name` - Can be aitoff, mollweide, arc, sinus, tan or mercator
  * @param {string} projection
  */
  setProjection(A) {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16), E = JA(A, h.__wbindgen_malloc, h.__wbindgen_realloc), o = MA;
      h.webclient_setProjection(I, this.__wbg_ptr, E, o);
      var g = f()[I / 4 + 0], C = f()[I / 4 + 1];
      if (C)
        throw AA(g);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Check whether the app is ready
  *
  * Aladin Lite is in a good state when the root tiles of the
  * HiPS chosen have all been retrieved and accessible for the GPU
  *
  * Surveys can be changed only if Aladin Lite is ready
  * @returns {boolean}
  */
  isReady() {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_isReady(I, this.__wbg_ptr);
      var A = f()[I / 4 + 0], g = f()[I / 4 + 1], C = f()[I / 4 + 2];
      if (C)
        throw AA(g);
      return A !== 0;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @returns {number}
  */
  getNOrder() {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_getNOrder(I, this.__wbg_ptr);
      var A = f()[I / 4 + 0], g = f()[I / 4 + 1], C = f()[I / 4 + 2];
      if (C)
        throw AA(g);
      return A;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Set new image surveys
  *
  * Send the image surveys to render inside the Aladin Lite view
  *
  * # Arguments
  *
  * * `surveys` - A list/array of survey. A survey is a javascript object
  * having the specific form. Please check the file in core/src/hips.rs to see
  * the different semantics accepted.
  *
  * # Examples
  *
  * ```javascript
  * let al = new Aladin.wasmLibs.core.WebClient(...);
  * const panstarrs = {
  *     properties: {
  *         url: "http://alasky.u-strasbg.fr/Pan-STARRS/DR1/r",
  *
  *         maxOrder: 11,
  *         frame: { label: "J2000", system: "J2000" },
  *         tileSize: 512,
  *         format: {
  *             FITSImage: {
  *                 bitpix: 16,
  *             }
  *         },
  *         minCutout: -0.15,
  *         maxCutout: 5,
  *     },
  *     color: {
  *         Grayscale2Colormap: {
  *             colormap: "RedTemperature",
  *             transfer: "asinh",
  *             reversed: false,
  *         }
  *     },
  * };
  * al.setImageSurveys([panstarrs]);
  * ```
  *
  * # Panics
  *
  * * If the surveys do not match SimpleHiPS type
  * * If the number of surveys is greater than 4. For the moment, due to the limitations
  *   of WebGL2 texture units on some architectures, the total number of surveys rendered is
  *   limited to 4.
  * @param {any} hips
  */
  addImageSurvey(A) {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_addImageSurvey(I, this.__wbg_ptr, x(A));
      var g = f()[I / 4 + 0], C = f()[I / 4 + 1];
      if (C)
        throw AA(g);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {any} fits_cfg
  * @returns {Promise<any>}
  */
  addImageFITS(A) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_addImageFITS(E, this.__wbg_ptr, x(A));
      var g = f()[E / 4 + 0], C = f()[E / 4 + 1], I = f()[E / 4 + 2];
      if (I)
        throw AA(C);
      return AA(g);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {string} layer
  */
  removeLayer(A) {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16), E = JA(A, h.__wbindgen_malloc, h.__wbindgen_realloc), o = MA;
      h.webclient_removeLayer(I, this.__wbg_ptr, E, o);
      var g = f()[I / 4 + 0], C = f()[I / 4 + 1];
      if (C)
        throw AA(g);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {string} layer
  * @param {string} new_layer
  */
  renameLayer(A, g) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16), o = JA(A, h.__wbindgen_malloc, h.__wbindgen_realloc), t = MA, s = JA(g, h.__wbindgen_malloc, h.__wbindgen_realloc), M = MA;
      h.webclient_renameLayer(E, this.__wbg_ptr, o, t, s, M);
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      if (I)
        throw AA(C);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {string} first_layer
  * @param {string} second_layer
  */
  swapLayers(A, g) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16), o = JA(A, h.__wbindgen_malloc, h.__wbindgen_realloc), t = MA, s = JA(g, h.__wbindgen_malloc, h.__wbindgen_realloc), M = MA;
      h.webclient_swapLayers(E, this.__wbg_ptr, o, t, s, M);
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      if (I)
        throw AA(C);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {string} past_url
  * @param {string} new_url
  */
  setHiPSUrl(A, g) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16), o = JA(A, h.__wbindgen_malloc, h.__wbindgen_realloc), t = MA, s = JA(g, h.__wbindgen_malloc, h.__wbindgen_realloc), M = MA;
      h.webclient_setHiPSUrl(E, this.__wbg_ptr, o, t, s, M);
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      if (I)
        throw AA(C);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {string} layer
  * @returns {ImageMetadata}
  */
  getImageMetadata(A) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16), o = JA(A, h.__wbindgen_malloc, h.__wbindgen_realloc), t = MA;
      h.webclient_getImageMetadata(E, this.__wbg_ptr, o, t);
      var g = f()[E / 4 + 0], C = f()[E / 4 + 1], I = f()[E / 4 + 2];
      if (I)
        throw AA(C);
      return vI.__wrap(g);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {string} layer
  * @param {any} meta
  */
  setImageMetadata(A, g) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16), o = JA(A, h.__wbindgen_malloc, h.__wbindgen_realloc), t = MA;
      h.webclient_setImageMetadata(E, this.__wbg_ptr, o, t, x(g));
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      if (I)
        throw AA(C);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {string} past_url
  * @param {string} new_url
  */
  setImageSurveyUrl(A, g) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16), o = JA(A, h.__wbindgen_malloc, h.__wbindgen_realloc), t = MA, s = JA(g, h.__wbindgen_malloc, h.__wbindgen_realloc), M = MA;
      h.webclient_setImageSurveyUrl(E, this.__wbg_ptr, o, t, s, M);
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      if (I)
        throw AA(C);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {any} color
  */
  setBackgroundColor(A) {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_setBackgroundColor(I, this.__wbg_ptr, x(A));
      var g = f()[I / 4 + 0], C = f()[I / 4 + 1];
      if (C)
        throw AA(g);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Set the equatorial grid color
  *
  * # Arguments
  *
  * * `red` - Red amount (between 0.0 and 1.0)
  * * `green` - Green amount (between 0.0 and 1.0)
  * * `blue` - Blue amount (between 0.0 and 1.0)
  * * `alpha` - Alpha amount (between 0.0 and 1.0)
  * @param {any} cfg
  */
  setGridConfig(A) {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_setGridConfig(I, this.__wbg_ptr, x(A));
      var g = f()[I / 4 + 0], C = f()[I / 4 + 1];
      if (C)
        throw AA(g);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Set the coordinate system for the view
  *
  * # Arguments
  *
  * * `coo_system` - The coordinate system
  * @param {number} coo_system
  */
  setCooSystem(A) {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_setCooSystem(I, this.__wbg_ptr, A);
      var g = f()[I / 4 + 0], C = f()[I / 4 + 1];
      if (C)
        throw AA(g);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Get the field of the view in degrees
  * @returns {number}
  */
  getFieldOfView() {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_getFieldOfView(I, this.__wbg_ptr);
      var A = AI()[I / 8 + 0], g = f()[I / 4 + 2], C = f()[I / 4 + 3];
      if (C)
        throw AA(g);
      return A;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Set the field of view
  *
  * # Arguments
  *
  * * `fov` - The field of view in degrees
  * @param {number} fov
  */
  setFieldOfView(A) {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_setFieldOfView(I, this.__wbg_ptr, A);
      var g = f()[I / 4 + 0], C = f()[I / 4 + 1];
      if (C)
        throw AA(g);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Set the absolute orientation of the view
  *
  * # Arguments
  *
  * * `theta` - The rotation angle in degrees
  * @param {number} theta
  */
  setRotationAroundCenter(A) {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_setRotationAroundCenter(I, this.__wbg_ptr, A);
      var g = f()[I / 4 + 0], C = f()[I / 4 + 1];
      if (C)
        throw AA(g);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Get the absolute orientation angle of the view
  * @returns {number}
  */
  getRotationAroundCenter() {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_getRotationAroundCenter(I, this.__wbg_ptr);
      var A = AI()[I / 8 + 0], g = f()[I / 4 + 2], C = f()[I / 4 + 3];
      if (C)
        throw AA(g);
      return A;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Get if the longitude axis is reversed
  * @returns {boolean}
  */
  getLongitudeReversed() {
    return h.webclient_getLongitudeReversed(this.__wbg_ptr) !== 0;
  }
  /**
  * Get the field of view angle value when the view is zoomed out to its maximum
  *
  * This method is dependent of the projection currently set.
  * All sky projections should return 360 degrees whereas
  * the sinus would be 180 degrees.
  * @returns {number}
  */
  getMaxFieldOfView() {
    return h.webclient_getMaxFieldOfView(this.__wbg_ptr);
  }
  /**
  * Get the clip zoom factor of the view
  *
  * This factor is deduced from the field of view angle.
  * It is a constant which when multiplied to the screen coordinates
  * gives the coordinates in clipping space.
  * @returns {number}
  */
  getClipZoomFactor() {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_getClipZoomFactor(I, this.__wbg_ptr);
      var A = AI()[I / 8 + 0], g = f()[I / 4 + 2], C = f()[I / 4 + 3];
      if (C)
        throw AA(g);
      return A;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Set the center of the view in ICRS coosys
  *
  * The core works in ICRS system so
  * the location must be given in this system
  *
  * # Arguments
  *
  * * `lon` - A longitude in degrees
  * * `lat` - A latitude in degrees
  * @param {number} lon
  * @param {number} lat
  */
  setCenter(A, g) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_setCenter(E, this.__wbg_ptr, A, g);
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      if (I)
        throw AA(C);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Get the center of the view
  *
  * This returns a javascript array of size 2.
  * The first component is the longitude, the second one is the latitude.
  * The angles are given in degrees.
  * @returns {Float64Array}
  */
  getCenter() {
    try {
      const o = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_getCenter(o, this.__wbg_ptr);
      var A = f()[o / 4 + 0], g = f()[o / 4 + 1], C = f()[o / 4 + 2], I = f()[o / 4 + 3];
      if (I)
        throw AA(C);
      var E = Pg(A, g).slice();
      return h.__wbindgen_free(A, g * 8), E;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Rest the north pole orientation to the top of the screen
  */
  resetNorthOrientation() {
    h.webclient_resetNorthOrientation(this.__wbg_ptr);
  }
  /**
  * Go from a location to another one
  *
  * # Arguments
  *
  * * `s1x` - The x screen coordinate in pixels of the starting point
  * * `s1y` - The y screen coordinate in pixels of the starting point
  * * `s2x` - The x screen coordinate in pixels of the goal point
  * * `s2y` - The y screen coordinate in pixels of the goal point
  * @param {number} s1x
  * @param {number} s1y
  * @param {number} s2x
  * @param {number} s2y
  */
  goFromTo(A, g, C, I) {
    try {
      const t = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_goFromTo(t, this.__wbg_ptr, A, g, C, I);
      var E = f()[t / 4 + 0], o = f()[t / 4 + 1];
      if (o)
        throw AA(E);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * View frame to ICRS/J2000 coosys conversion
  *
  * Coordinates must be given in the ICRS coo system
  *
  * # Arguments
  *
  * * `lon` - A longitude in degrees
  * * `lat` - A latitude in degrees
  * @param {number} lon
  * @param {number} lat
  * @returns {Float64Array}
  */
  viewToICRSCooSys(A, g) {
    try {
      const o = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_viewToICRSCooSys(o, this.__wbg_ptr, A, g);
      var C = f()[o / 4 + 0], I = f()[o / 4 + 1], E = Pg(C, I).slice();
      return h.__wbindgen_free(C, I * 8), E;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * World to screen projection
  *
  * Coordinates must be given in the ICRS coo system
  *
  * # Arguments
  *
  * * `lon` - A longitude in degrees
  * * `lat` - A latitude in degrees
  * @param {number} lon
  * @param {number} lat
  * @returns {Float64Array | undefined}
  */
  worldToScreen(A, g) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_worldToScreen(E, this.__wbg_ptr, A, g);
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      let o;
      return C !== 0 && (o = Pg(C, I).slice(), h.__wbindgen_free(C, I * 8)), o;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {number} x
  * @param {number} y
  * @returns {Float64Array}
  */
  screenToClip(A, g) {
    try {
      const o = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_screenToClip(o, this.__wbg_ptr, A, g);
      var C = f()[o / 4 + 0], I = f()[o / 4 + 1], E = Pg(C, I).slice();
      return h.__wbindgen_free(C, I * 8), E;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {Float64Array} lon
  * @param {Float64Array} lat
  * @returns {Float64Array}
  */
  worldToScreenVec(A, g) {
    try {
      const o = h.__wbindgen_add_to_stack_pointer(-16), t = DQ(A, h.__wbindgen_malloc), s = MA, M = DQ(g, h.__wbindgen_malloc), N = MA;
      h.webclient_worldToScreenVec(o, this.__wbg_ptr, t, s, M, N);
      var C = f()[o / 4 + 0], I = f()[o / 4 + 1], E = Pg(C, I).slice();
      return h.__wbindgen_free(C, I * 8), E;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {any} _catalog
  */
  setCatalog(A) {
    try {
      h.webclient_setCatalog(this.__wbg_ptr, aQ(A));
    } finally {
      Mg[NI++] = void 0;
    }
  }
  /**
  * Screen to world unprojection
  *
  * # Arguments
  *
  * * `pos_x` - The x screen coordinate in pixels
  * * `pos_y` - The y screen coordinate in pixels
  * @param {number} pos_x
  * @param {number} pos_y
  * @returns {Float64Array | undefined}
  */
  screenToWorld(A, g) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_screenToWorld(E, this.__wbg_ptr, A, g);
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      let o;
      return C !== 0 && (o = Pg(C, I).slice(), h.__wbindgen_free(C, I * 8)), o;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Signal the backend when the left mouse button has been released.
  *
  * This is useful for beginning inerting.
  * @param {number} sx
  * @param {number} sy
  */
  releaseLeftButtonMouse(A, g) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_releaseLeftButtonMouse(E, this.__wbg_ptr, A, g);
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      if (I)
        throw AA(C);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Signal the backend when the left mouse button has been pressed.
  * @param {number} sx
  * @param {number} sy
  */
  pressLeftMouseButton(A, g) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_pressLeftMouseButton(E, this.__wbg_ptr, A, g);
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      if (I)
        throw AA(C);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Signal the backend when the left mouse button has been pressed.
  * @param {number} s1x
  * @param {number} s1y
  * @param {number} s2x
  * @param {number} s2y
  */
  moveMouse(A, g, C, I) {
    try {
      const t = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_moveMouse(t, this.__wbg_ptr, A, g, C, I);
      var E = f()[t / 4 + 0], o = f()[t / 4 + 1];
      if (o)
        throw AA(E);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Add a catalog rendered as a heatmap.
  *
  * # Arguments
  *
  * * `name_catalog` - The name of the catalog
  * * `data` - The list of the catalog sources.
  * * `colormap` - The name of the colormap. Check out the list of possible colormaps names `getAvailableColormapList`.
  * @param {string} name_catalog
  * @param {any} data
  * @param {string} colormap
  */
  addCatalog(A, g, C) {
    try {
      const o = h.__wbindgen_add_to_stack_pointer(-16), t = JA(A, h.__wbindgen_malloc, h.__wbindgen_realloc), s = MA, M = JA(C, h.__wbindgen_malloc, h.__wbindgen_realloc), N = MA;
      h.webclient_addCatalog(o, this.__wbg_ptr, t, s, x(g), M, N);
      var I = f()[o / 4 + 0], E = f()[o / 4 + 1];
      if (E)
        throw AA(I);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Set the catalog heatmap colormap
  *
  * # Arguments
  *
  * * `name_catalog` - The name of the catalog to apply this change to
  * * `colormap` - The name of the colormap. Check out the list of possible colormaps names `getAvailableColormapList`.
  *
  * # Panics
  *
  * If the catalog has not been found
  * @returns {boolean}
  */
  isCatalogLoaded() {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_isCatalogLoaded(I, this.__wbg_ptr);
      var A = f()[I / 4 + 0], g = f()[I / 4 + 1], C = f()[I / 4 + 2];
      if (C)
        throw AA(g);
      return A !== 0;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Set the catalog heatmap opacity
  *
  * # Arguments
  *
  * * `name_catalog` - The name of the catalog to apply this change to
  * * `opacity` - The opacity factor (between 0.0 and 1.0)
  *
  * # Panics
  *
  * If the catalog has not been found
  * @param {string} name_catalog
  * @param {number} opacity
  */
  setCatalogOpacity(A, g) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16), o = JA(A, h.__wbindgen_malloc, h.__wbindgen_realloc), t = MA;
      h.webclient_setCatalogOpacity(E, this.__wbg_ptr, o, t, g);
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      if (I)
        throw AA(C);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Set the kernel strength for the catalog heatmap rendering
  *
  * # Arguments
  *
  * * `name_catalog` - The name of the catalog to apply this change to
  * * `strength` - The strength of the kernel
  *
  * # Panics
  *
  * If the catalog has not been found
  * @param {string} name_catalog
  * @param {number} strength
  */
  setCatalogKernelStrength(A, g) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16), o = JA(A, h.__wbindgen_malloc, h.__wbindgen_realloc), t = MA;
      h.webclient_setCatalogKernelStrength(E, this.__wbg_ptr, o, t, g);
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      if (I)
        throw AA(C);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Project a line to the screen
  *
  * # Returns
  *
  * A list of xy screen coordinates defining the projected line.
  * The algorithm involved is recursive and can return different number of
  * control points depending on the projection used and therefore
  * the deformation of the line.
  *
  * # Arguments
  *
  * * `lon1` - The longitude in degrees of the starting line point
  * * `lat1` - The latitude in degrees of the starting line point
  * * `lon2` - The longitude in degrees of the ending line point
  * * `lat2` - The latitude in degrees of the ending line point
  * Get the list of colormap supported
  *
  * This list must be updated whenever a new colormap is added
  * in core/img/colormaps/colormaps.png
  * @returns {any[]}
  */
  getAvailableColormapList() {
    try {
      const o = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_getAvailableColormapList(o, this.__wbg_ptr);
      var A = f()[o / 4 + 0], g = f()[o / 4 + 1], C = f()[o / 4 + 2], I = f()[o / 4 + 3];
      if (I)
        throw AA(C);
      var E = ZE(A, g).slice();
      return h.__wbindgen_free(A, g * 4), E;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {string} label
  * @param {any[]} hex_colors
  */
  createCustomColormap(A, g) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16), o = JA(A, h.__wbindgen_malloc, h.__wbindgen_realloc), t = MA, s = TE(g, h.__wbindgen_malloc), M = MA;
      h.webclient_createCustomColormap(E, this.__wbg_ptr, o, t, s, M);
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      if (I)
        throw AA(C);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Get the image canvas where the webgl rendering is done
  * @returns {object | undefined}
  */
  canvas() {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_canvas(I, this.__wbg_ptr);
      var A = f()[I / 4 + 0], g = f()[I / 4 + 1], C = f()[I / 4 + 2];
      if (C)
        throw AA(g);
      return AA(A);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * Read the pixel value
  *
  * The current implementation only returns the pixel value
  * of the first survey of the `layer` specified.
  *
  * # Returns
  *
  * - An array of 3 items (rgb) for JPG tiles
  * - An array of 4 items (rgba) for PNG tiles
  * - A single value for FITS tiles
  *
  * # Arguments
  *
  * * `x` - The x screen coordinate in pixels
  * * `y` - The y screen coordinate in pixels
  * * `base_url` - The base url of the survey identifying it
  * @param {number} x
  * @param {number} y
  * @param {string} layer
  * @returns {any}
  */
  readPixel(A, g, C) {
    try {
      const t = h.__wbindgen_add_to_stack_pointer(-16), s = JA(C, h.__wbindgen_malloc, h.__wbindgen_realloc), M = MA;
      h.webclient_readPixel(t, this.__wbg_ptr, A, g, s, M);
      var I = f()[t / 4 + 0], E = f()[t / 4 + 1], o = f()[t / 4 + 2];
      if (o)
        throw AA(E);
      return AA(I);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {number} depth
  * @returns {any}
  */
  getVisibleCells(A) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16);
      h.webclient_getVisibleCells(E, this.__wbg_ptr, A);
      var g = f()[E / 4 + 0], C = f()[E / 4 + 1], I = f()[E / 4 + 2];
      if (I)
        throw AA(C);
      return AA(g);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @returns {boolean}
  */
  isRendering() {
    return h.webclient_isRendering(this.__wbg_ptr) !== 0;
  }
  /**
  * @param {MOC} params
  * @param {any} data
  */
  addJSONMoc(A, g) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16);
      kg(A, ng), h.webclient_addJSONMoc(E, this.__wbg_ptr, A.__wbg_ptr, aQ(g));
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      if (I)
        throw AA(C);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16), Mg[NI++] = void 0;
    }
  }
  /**
  * @param {string} s
  * @returns {any}
  */
  parseVOTable(A) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16), o = JA(A, h.__wbindgen_malloc, h.__wbindgen_realloc), t = MA;
      h.webclient_parseVOTable(E, this.__wbg_ptr, o, t);
      var g = f()[E / 4 + 0], C = f()[E / 4 + 1], I = f()[E / 4 + 2];
      if (I)
        throw AA(C);
      return AA(g);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {MOC} params
  * @param {Uint8Array} data
  */
  addFITSMoc(A, g) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16);
      kg(A, ng);
      const o = lQ(g, h.__wbindgen_malloc), t = MA;
      h.webclient_addFITSMoc(E, this.__wbg_ptr, A.__wbg_ptr, o, t);
      var C = f()[E / 4 + 0], I = f()[E / 4 + 1];
      if (I)
        throw AA(C);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {MOC} params
  */
  removeMoc(A) {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16);
      kg(A, ng), h.webclient_removeMoc(I, this.__wbg_ptr, A.__wbg_ptr);
      var g = f()[I / 4 + 0], C = f()[I / 4 + 1];
      if (C)
        throw AA(g);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {MOC} cfg
  */
  setMocParams(A) {
    try {
      const I = h.__wbindgen_add_to_stack_pointer(-16);
      kg(A, ng), h.webclient_setMocParams(I, this.__wbg_ptr, A.__wbg_ptr);
      var g = f()[I / 4 + 0], C = f()[I / 4 + 1];
      if (C)
        throw AA(g);
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {MOC} _params
  * @param {number} _lon
  * @param {number} _lat
  * @returns {boolean}
  */
  mocContains(A, g, C) {
    try {
      const t = h.__wbindgen_add_to_stack_pointer(-16);
      kg(A, ng), h.webclient_mocContains(t, this.__wbg_ptr, A.__wbg_ptr, g, C);
      var I = f()[t / 4 + 0], E = f()[t / 4 + 1], o = f()[t / 4 + 2];
      if (o)
        throw AA(E);
      return I !== 0;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
  /**
  * @param {MOC} _params
  * @returns {number}
  */
  mocSkyFraction(A) {
    try {
      const E = h.__wbindgen_add_to_stack_pointer(-16);
      kg(A, ng), h.webclient_mocSkyFraction(E, this.__wbg_ptr, A.__wbg_ptr);
      var g = kI()[E / 4 + 0], C = f()[E / 4 + 1], I = f()[E / 4 + 2];
      if (I)
        throw AA(C);
      return g;
    } finally {
      h.__wbindgen_add_to_stack_pointer(16);
    }
  }
}
async function oi(B, A) {
  if (typeof Response == "function" && B instanceof Response) {
    if (typeof WebAssembly.instantiateStreaming == "function")
      try {
        return await WebAssembly.instantiateStreaming(B, A);
      } catch (C) {
        if (B.headers.get("Content-Type") != "application/wasm")
          console.warn("`WebAssembly.instantiateStreaming` failed because your server does not serve wasm with `application/wasm` MIME type. Falling back to `WebAssembly.instantiate` which is slower. Original error:\n", C);
        else
          throw C;
      }
    const g = await B.arrayBuffer();
    return await WebAssembly.instantiate(g, A);
  } else {
    const g = await WebAssembly.instantiate(B, A);
    return g instanceof WebAssembly.Instance ? { instance: g, module: B } : g;
  }
}
function JQ() {
  const B = {};
  return B.wbg = {}, B.wbg.__wbg_hexToRgba_47dea88f7840631a = function(A, g) {
    let C, I;
    try {
      C = A, I = g;
      const E = _A.hexToRgba(dA(A, g));
      return x(E);
    } finally {
      h.__wbindgen_free(C, I, 1);
    }
  }, B.wbg.__wbg_getwithrefkey_15c62c2b8546208d = function(A, g) {
    const C = r(A)[r(g)];
    return x(C);
  }, B.wbg.__wbg_warn_d60e832f9882c1b2 = function(A) {
    console.warn(r(A));
  }, B.wbg.__wbindgen_string_new = function(A, g) {
    const C = dA(A, g);
    return x(C);
  }, B.wbg.__wbg_length_fff51ee6522a1a18 = function(A) {
    return r(A).length;
  }, B.wbg.__wbg_getElementById_cc0e0d931b0d9a28 = function(A, g, C) {
    const I = r(A).getElementById(dA(g, C));
    return FA(I) ? 0 : x(I);
  }, B.wbg.__wbg_getElementsByClassName_e29c7f3f729baad5 = function(A, g, C) {
    const I = r(A).getElementsByClassName(dA(g, C));
    return x(I);
  }, B.wbg.__wbg_parse_670c19d4e984792e = function() {
    return kA(function(A, g) {
      const C = JSON.parse(dA(A, g));
      return x(C);
    }, arguments);
  }, B.wbg.__wbg_getContext_6d1f155bb5c1096a = function() {
    return kA(function(A, g, C, I) {
      const E = r(A).getContext(dA(g, C), r(I));
      return FA(E) ? 0 : x(E);
    }, arguments);
  }, B.wbg.__wbg_instanceof_WebGl2RenderingContext_f921526c513bf717 = function(A) {
    let g;
    try {
      g = r(A) instanceof WebGL2RenderingContext;
    } catch {
      g = !1;
    }
    return g;
  }, B.wbg.__wbg_getExtension_77909f6d51d49d4d = function() {
    return kA(function(A, g, C) {
      const I = r(A).getExtension(dA(g, C));
      return FA(I) ? 0 : x(I);
    }, arguments);
  }, B.wbg.__wbg_newnoargs_581967eacc0e2604 = function(A, g) {
    const C = new Function(dA(A, g));
    return x(C);
  }, B.wbg.__wbg_style_3801009b2339aa94 = function(A) {
    const g = r(A).style;
    return x(g);
  }, B.wbg.__wbg_setwidth_a667a942dba6656e = function(A, g) {
    r(A).width = g >>> 0;
  }, B.wbg.__wbg_setheight_a747d440760fe5aa = function(A, g) {
    r(A).height = g >>> 0;
  }, B.wbg.__wbg_getwithrefkey_5e6d9547403deab8 = function(A, g) {
    const C = r(A)[r(g)];
    return x(C);
  }, B.wbg.__wbg_get_44be0491f933a435 = function(A, g) {
    const C = r(A)[g >>> 0];
    return x(C);
  }, B.wbg.__wbg_new_d258248ed531ff54 = function(A, g) {
    const C = new Error(dA(A, g));
    return x(C);
  }, B.wbg.__wbg_new_43f1b47c28813cbd = function(A, g) {
    try {
      var C = { a: A, b: g }, I = (o, t) => {
        const s = C.a;
        C.a = 0;
        try {
          return WE(s, C.b, o, t);
        } finally {
          C.a = s;
        }
      };
      const E = new Promise(I);
      return x(E);
    } finally {
      C.a = C.b = 0;
    }
  }, B.wbg.__wbg_finally_d33c80d30b8ede4e = function(A, g) {
    const C = r(A).finally(r(g));
    return x(C);
  }, B.wbg.__wbindgen_object_clone_ref = function(A) {
    const g = r(A);
    return x(g);
  }, B.wbg.__wbg_rgbToHex_c6eef179c7d74290 = function(A, g, C, I) {
    const E = _A.rgbToHex(g, C, I), o = JA(E, h.__wbindgen_malloc, h.__wbindgen_realloc), t = MA;
    f()[A / 4 + 1] = t, f()[A / 4 + 0] = o;
  }, B.wbg.__wbindgen_object_drop_ref = function(A) {
    AA(A);
  }, B.wbg.__wbg_checkFramebufferStatus_d5eaeabeb7ea0712 = function(A, g) {
    return r(A).checkFramebufferStatus(g >>> 0);
  }, B.wbg.__wbg_newwithlength_e5d69174d6984cd7 = function(A) {
    const g = new Uint8Array(A >>> 0);
    return x(g);
  }, B.wbg.__wbg_new_898a68150f225f2e = function() {
    const A = new Array();
    return x(A);
  }, B.wbg.__wbg_newwithlength_68d29ab115d0099c = function(A) {
    const g = new Float32Array(A >>> 0);
    return x(g);
  }, B.wbg.__wbg_length_d7327c75a759af37 = function(A) {
    return r(A).length;
  }, B.wbg.__wbindgen_memory = function() {
    const A = h.memory;
    return x(A);
  }, B.wbg.__wbg_buffer_085ec1f694018c4f = function(A) {
    const g = r(A).buffer;
    return x(g);
  }, B.wbg.__wbg_new_d086a66d1c264b3f = function(A) {
    const g = new Float32Array(r(A));
    return x(g);
  }, B.wbg.__wbg_set_6146c51d49a2c0df = function(A, g, C) {
    r(A).set(r(g), C >>> 0);
  }, B.wbg.__wbg_newwithlength_63d69a995741f1a4 = function(A) {
    const g = new Int16Array(A >>> 0);
    return x(g);
  }, B.wbg.__wbg_length_86d47c4a8fd89d6f = function(A) {
    return r(A).length;
  }, B.wbg.__wbg_new_02e18b6397a0f2d8 = function(A) {
    const g = new Int16Array(r(A));
    return x(g);
  }, B.wbg.__wbg_set_96d7d44309334124 = function(A, g, C) {
    r(A).set(r(g), C >>> 0);
  }, B.wbg.__wbindgen_number_new = function(A) {
    return x(A);
  }, B.wbg.__wbg_newwithlength_17ac8c88a224411a = function(A) {
    const g = new Int32Array(A >>> 0);
    return x(g);
  }, B.wbg.__wbg_length_10541883ff7175cc = function(A) {
    return r(A).length;
  }, B.wbg.__wbg_new_a0af68041688e8fd = function(A) {
    const g = new Int32Array(r(A));
    return x(g);
  }, B.wbg.__wbg_set_8fe6d6fc52f887cb = function(A, g, C) {
    r(A).set(r(g), C >>> 0);
  }, B.wbg.__wbg_width_2931aaedd21f1fff = function(A) {
    return r(A).width;
  }, B.wbg.__wbg_height_0d36fbbeb60b0661 = function(A) {
    return r(A).height;
  }, B.wbg.__wbindgen_number_get = function(A, g) {
    const C = r(g), I = typeof C == "number" ? C : void 0;
    AI()[A / 8 + 1] = FA(I) ? 0 : I, f()[A / 4 + 0] = !FA(I);
  }, B.wbg.__wbg_new_b51585de1b234aff = function() {
    const A = new Object();
    return x(A);
  }, B.wbg.__wbg_set_502d29070ea18557 = function(A, g, C) {
    r(A)[g >>> 0] = AA(C);
  }, B.wbg.__wbg_stringify_e25465938f3f611f = function() {
    return kA(function(A) {
      const g = JSON.stringify(r(A));
      return x(g);
    }, arguments);
  }, B.wbg.__wbindgen_string_get = function(A, g) {
    const C = r(g), I = typeof C == "string" ? C : void 0;
    var E = FA(I) ? 0 : JA(I, h.__wbindgen_malloc, h.__wbindgen_realloc), o = MA;
    f()[A / 4 + 1] = o, f()[A / 4 + 0] = E;
  }, B.wbg.__wbg_set_bedc3d02d0f05eb0 = function(A, g, C) {
    const I = r(A).set(r(g), r(C));
    return x(I);
  }, B.wbg.__wbg_set_841ac57cff3d672b = function(A, g, C) {
    r(A)[AA(g)] = AA(C);
  }, B.wbg.__wbindgen_is_string = function(A) {
    return typeof r(A) == "string";
  }, B.wbg.__wbindgen_is_object = function(A) {
    const g = r(A);
    return typeof g == "object" && g !== null;
  }, B.wbg.__wbg_entries_e51f29c7bba0c054 = function(A) {
    const g = Object.entries(r(A));
    return x(g);
  }, B.wbg.__wbindgen_boolean_get = function(A) {
    const g = r(A);
    return typeof g == "boolean" ? g ? 1 : 0 : 2;
  }, B.wbg.__wbindgen_is_undefined = function(A) {
    return r(A) === void 0;
  }, B.wbg.__wbindgen_in = function(A, g) {
    return r(A) in r(g);
  }, B.wbg.__wbg_hexToRgb_196aacecf6bbc8ac = function(A, g) {
    let C, I;
    try {
      C = A, I = g;
      const E = _A.hexToRgb(dA(A, g));
      return x(E);
    } finally {
      h.__wbindgen_free(C, I, 1);
    }
  }, B.wbg.__wbg_createVertexArray_51d51e1e1e13e9f6 = function(A) {
    const g = r(A).createVertexArray();
    return FA(g) ? 0 : x(g);
  }, B.wbg.__wbg_deleteVertexArray_3e4f2e2ff7f05a19 = function(A, g) {
    r(A).deleteVertexArray(r(g));
  }, B.wbg.__wbg_drawElements_a9529eefaf2008bd = function(A, g, C, I, E) {
    r(A).drawElements(g >>> 0, C, I >>> 0, E);
  }, B.wbg.__wbg_vertexAttribPointer_316ffe2f0458fde7 = function(A, g, C, I, E, o, t) {
    r(A).vertexAttribPointer(g >>> 0, C, I >>> 0, E !== 0, o, t);
  }, B.wbg.__wbg_enableVertexAttribArray_8804480c2ea0bb72 = function(A, g) {
    r(A).enableVertexAttribArray(g >>> 0);
  }, B.wbg.__wbg_texSubImage2D_d2841ded12a8aa66 = function() {
    return kA(function(A, g, C, I, E, o, t, s, M, N) {
      r(A).texSubImage2D(g >>> 0, C, I, E, o, t, s >>> 0, M >>> 0, r(N));
    }, arguments);
  }, B.wbg.__wbg_deleteTexture_4fcfea73cd8f6214 = function(A, g) {
    r(A).deleteTexture(r(g));
  }, B.wbg.__wbg_texSubImage2D_b34eb18892e3f6e5 = function() {
    return kA(function(A, g, C, I, E, o, t, s) {
      r(A).texSubImage2D(g >>> 0, C, I, E, o >>> 0, t >>> 0, r(s));
    }, arguments);
  }, B.wbg.__wbg_texSubImage2D_ed22054fe33dee5e = function() {
    return kA(function(A, g, C, I, E, o, t, s) {
      r(A).texSubImage2D(g >>> 0, C, I, E, o >>> 0, t >>> 0, r(s));
    }, arguments);
  }, B.wbg.__wbg_new_8125e318e6245eed = function(A) {
    const g = new Uint8Array(r(A));
    return x(g);
  }, B.wbg.__wbg_texSubImage2D_2268f0b57770632e = function() {
    return kA(function(A, g, C, I, E, o, t, s) {
      r(A).texSubImage2D(g >>> 0, C, I, E, o >>> 0, t >>> 0, r(s));
    }, arguments);
  }, B.wbg.__wbg_createProgram_4eaf3b97b5747a62 = function(A) {
    const g = r(A).createProgram();
    return FA(g) ? 0 : x(g);
  }, B.wbg.__wbg_linkProgram_33998194075d71fb = function(A, g) {
    r(A).linkProgram(r(g));
  }, B.wbg.__wbg_getProgramInfoLog_b81bc53188e286fa = function(A, g, C) {
    const I = r(g).getProgramInfoLog(r(C));
    var E = FA(I) ? 0 : JA(I, h.__wbindgen_malloc, h.__wbindgen_realloc), o = MA;
    f()[A / 4 + 1] = o, f()[A / 4 + 0] = E;
  }, B.wbg.__wbg_getActiveUniform_78367ddc7339640b = function(A, g, C) {
    const I = r(A).getActiveUniform(r(g), C >>> 0);
    return FA(I) ? 0 : x(I);
  }, B.wbg.__wbg_name_ebae3a7e89367611 = function(A, g) {
    const C = r(g).name, I = JA(C, h.__wbindgen_malloc, h.__wbindgen_realloc), E = MA;
    f()[A / 4 + 1] = E, f()[A / 4 + 0] = I;
  }, B.wbg.__wbg_getUniformLocation_9f6eb60c560a347b = function(A, g, C, I) {
    const E = r(A).getUniformLocation(r(g), dA(C, I));
    return FA(E) ? 0 : x(E);
  }, B.wbg.__wbg_attachShader_47256b6b3d42a22e = function(A, g, C) {
    r(A).attachShader(r(g), r(C));
  }, B.wbg.__wbg_createShader_429776c9dd6fb87b = function(A, g) {
    const C = r(A).createShader(g >>> 0);
    return FA(C) ? 0 : x(C);
  }, B.wbg.__wbg_shaderSource_1cb7c64dc7d1a500 = function(A, g, C, I) {
    r(A).shaderSource(r(g), dA(C, I));
  }, B.wbg.__wbg_compileShader_6bf78b425d5c98e1 = function(A, g) {
    r(A).compileShader(r(g));
  }, B.wbg.__wbg_getShaderParameter_ac2727ae4fe7648e = function(A, g, C) {
    const I = r(A).getShaderParameter(r(g), C >>> 0);
    return x(I);
  }, B.wbg.__wbg_getShaderInfoLog_968b93e75477d725 = function(A, g, C) {
    const I = r(g).getShaderInfoLog(r(C));
    var E = FA(I) ? 0 : JA(I, h.__wbindgen_malloc, h.__wbindgen_realloc), o = MA;
    f()[A / 4 + 1] = o, f()[A / 4 + 0] = E;
  }, B.wbg.__wbg_activeTexture_799bf1387e911c27 = function(A, g) {
    r(A).activeTexture(g >>> 0);
  }, B.wbg.__wbg_cancel_7f202496da02cd45 = function(A) {
    const g = r(A).cancel();
    return x(g);
  }, B.wbg.__wbg_catch_64e0c7dcea0da34e = function(A, g) {
    const C = r(A).catch(r(g));
    return x(C);
  }, B.wbg.__wbg_body_b86f372950de5b7d = function(A) {
    const g = r(A).body;
    return FA(g) ? 0 : x(g);
  }, B.wbg.__wbg_instanceof_ReadableStream_723f1212419028dc = function(A) {
    let g;
    try {
      g = r(A) instanceof ReadableStream;
    } catch {
      g = !1;
    }
    return g;
  }, B.wbg.__wbg_getReader_dbfa3af83c1c7ed0 = function() {
    return kA(function(A, g) {
      const C = r(A).getReader(bI.__wrap(g));
      return x(C);
    }, arguments);
  }, B.wbg.__wbg_getReader_8ecba87d8003e950 = function() {
    return kA(function(A) {
      const g = r(A).getReader();
      return x(g);
    }, arguments);
  }, B.wbg.__wbg_getParameter_55b36a787dbbfb74 = function() {
    return kA(function(A, g) {
      const C = r(A).getParameter(g >>> 0);
      return x(C);
    }, arguments);
  }, B.wbg.__wbg_texImage2D_8dca6270a1399d43 = function() {
    return kA(function(A, g, C, I, E, o, t) {
      r(A).texImage2D(g >>> 0, C, I, E >>> 0, o >>> 0, r(t));
    }, arguments);
  }, B.wbg.__wbg_scissor_e8e41e1c0a9817c8 = function(A, g, C, I, E) {
    r(A).scissor(g, C, I, E);
  }, B.wbg.__wbg_instanceof_Response_fc4327dbfcdf5ced = function(A) {
    let g;
    try {
      g = r(A) instanceof Response;
    } catch {
      g = !1;
    }
    return g;
  }, B.wbg.__wbg_length_72e2208bbc0efc61 = function(A) {
    return r(A).length;
  }, B.wbg.__wbg_setonload_b4f5d9b15b0ee9d3 = function(A, g) {
    r(A).onload = r(g);
  }, B.wbg.__wbg_setonerror_acddd28c276005c1 = function(A, g) {
    r(A).onerror = r(g);
  }, B.wbg.__wbg_setsrc_fac5b9516fc69301 = function(A, g, C) {
    r(A).src = dA(g, C);
  }, B.wbg.__wbg_createElement_4891554b28d3388b = function() {
    return kA(function(A, g, C) {
      const I = r(A).createElement(dA(g, C));
      return x(I);
    }, arguments);
  }, B.wbg.__wbg_drawImage_9758fa4c03ab8fc8 = function() {
    return kA(function(A, g, C, I) {
      r(A).drawImage(r(g), C, I);
    }, arguments);
  }, B.wbg.__wbg_getImageData_cacdc8ba8433e1ff = function() {
    return kA(function(A, g, C, I, E) {
      const o = r(A).getImageData(g, C, I, E);
      return x(o);
    }, arguments);
  }, B.wbg.__wbg_data_eaf4962120932fdc = function(A, g) {
    const C = r(g).data, I = lQ(C, h.__wbindgen_malloc), E = MA;
    f()[A / 4 + 1] = E, f()[A / 4 + 0] = I;
  }, B.wbg.__wbg_ok_e3d8d84e630fd064 = function(A) {
    return r(A).ok;
  }, B.wbg.__wbg_byteLength_47d11fa79875dee3 = function(A) {
    return r(A).byteLength;
  }, B.wbg.__wbg_read_88c96573fc8b3b01 = function(A) {
    const g = r(A).read();
    return x(g);
  }, B.wbg.__wbg_done_76252d32deca186b = function(A) {
    return r(A).done;
  }, B.wbg.__wbg_value_ff3741eb46856618 = function(A) {
    const g = r(A).value;
    return x(g);
  }, B.wbg.__wbg_subarray_13db269f57aa838d = function(A, g, C) {
    const I = r(A).subarray(g >>> 0, C >>> 0);
    return x(I);
  }, B.wbg.__wbg_read_506484ad103ddec0 = function(A, g) {
    const C = r(A).read(r(g));
    return x(C);
  }, B.wbg.__wbg_done_d22551bb874c751e = function(A) {
    return r(A).done;
  }, B.wbg.__wbg_toString_a8e343996af880e9 = function(A) {
    const g = r(A).toString();
    return x(g);
  }, B.wbg.__wbg_value_001280e0aafa4bda = function(A) {
    const g = r(A).value;
    return FA(g) ? 0 : x(g);
  }, B.wbg.__wbg_buffer_f5b7059c439f330d = function(A) {
    const g = r(A).buffer;
    return x(g);
  }, B.wbg.__wbg_innerWidth_ebe07ce5463ff293 = function() {
    return kA(function(A) {
      const g = r(A).innerWidth;
      return x(g);
    }, arguments);
  }, B.wbg.__wbg_innerHeight_2dd06d8cf68f1d7d = function() {
    return kA(function(A) {
      const g = r(A).innerHeight;
      return x(g);
    }, arguments);
  }, B.wbg.__wbg_devicePixelRatio_f9de7bddca0eaf20 = function(A) {
    return r(A).devicePixelRatio;
  }, B.wbg.__wbg_width_e64ae54f1609bb76 = function(A) {
    return r(A).width;
  }, B.wbg.__wbg_height_5ee3e7570341fe45 = function(A) {
    return r(A).height;
  }, B.wbg.__wbg_getElementsByClassName_173ed29f6818e275 = function(A, g, C) {
    const I = r(A).getElementsByClassName(dA(g, C));
    return x(I);
  }, B.wbg.__wbg_clearColor_7a7d04702f7e38e5 = function(A, g, C, I, E) {
    r(A).clearColor(g, C, I, E);
  }, B.wbg.__wbg_clear_2db2efe323bfdf68 = function(A, g) {
    r(A).clear(g >>> 0);
  }, B.wbg.__wbg_uniform3f_8364a0959b6c1570 = function(A, g, C, I, E) {
    r(A).uniform3f(r(g), C, I, E);
  }, B.wbg.__wbg_uniform4f_a9fd337d4b07f595 = function(A, g, C, I, E, o) {
    r(A).uniform4f(r(g), C, I, E, o);
  }, B.wbg.__wbg_from_d7c216d4616bb368 = function(A) {
    const g = Array.from(r(A));
    return x(g);
  }, B.wbg.__wbg_performance_2c295061c8b01e0b = function(A) {
    const g = r(A).performance;
    return FA(g) ? 0 : x(g);
  }, B.wbg.__wbg_now_0cfdc90c97d0c24b = function(A) {
    return r(A).now();
  }, B.wbg.__wbg_subarray_6814da0bb52e50aa = function(A, g, C) {
    const I = r(A).subarray(g >>> 0, C >>> 0);
    return x(I);
  }, B.wbg.__wbg_instanceof_CanvasRenderingContext2d_bc0a6635c96eca9b = function(A) {
    let g;
    try {
      g = r(A) instanceof CanvasRenderingContext2D;
    } catch {
      g = !1;
    }
    return g;
  }, B.wbg.__wbindgen_error_new = function(A, g) {
    const C = new Error(dA(A, g));
    return x(C);
  }, B.wbg.__wbg_uniformMatrix4fv_766b5ba343983038 = function(A, g, C, I, E) {
    r(A).uniformMatrix4fv(r(g), C !== 0, jE(I, E));
  }, B.wbg.__wbg_clearRect_517d3360d8be8a55 = function(A, g, C, I, E) {
    r(A).clearRect(g, C, I, E);
  }, B.wbg.__wbg_setfont_3d2943420717ac87 = function(A, g, C) {
    r(A).font = dA(g, C);
  }, B.wbg.__wbg_setfillStyle_401fa583a1c8863c = function(A, g) {
    r(A).fillStyle = r(g);
  }, B.wbg.__wbg_save_cdcca9591f027e80 = function(A) {
    r(A).save();
  }, B.wbg.__wbg_translate_6156c415b2341ec9 = function() {
    return kA(function(A, g, C) {
      r(A).translate(g, C);
    }, arguments);
  }, B.wbg.__wbg_rotate_5e5e156d10e67adc = function() {
    return kA(function(A, g) {
      r(A).rotate(g);
    }, arguments);
  }, B.wbg.__wbg_settextAlign_2601a967e1935930 = function(A, g, C) {
    r(A).textAlign = dA(g, C);
  }, B.wbg.__wbg_fillText_ba4313e6835ce7ea = function() {
    return kA(function(A, g, C, I, E) {
      r(A).fillText(dA(g, C), I, E);
    }, arguments);
  }, B.wbg.__wbg_restore_890c3582852dbadf = function(A) {
    r(A).restore();
  }, B.wbg.__wbg_new_abda76e883ba8a5f = function() {
    const A = new Error();
    return x(A);
  }, B.wbg.__wbg_stack_658279fe44541cf6 = function(A, g) {
    const C = r(g).stack, I = JA(C, h.__wbindgen_malloc, h.__wbindgen_realloc), E = MA;
    f()[A / 4 + 1] = E, f()[A / 4 + 0] = I;
  }, B.wbg.__wbg_error_f851667af71bcfc6 = function(A, g) {
    let C, I;
    try {
      C = A, I = g, console.error(dA(A, g));
    } finally {
      h.__wbindgen_free(C, I, 1);
    }
  }, B.wbg.__wbg_next_ddb3312ca1c4e32a = function() {
    return kA(function(A) {
      const g = r(A).next();
      return x(g);
    }, arguments);
  }, B.wbg.__wbg_done_5c1f01fb660d73b5 = function(A) {
    return r(A).done;
  }, B.wbg.__wbg_value_1695675138684bd5 = function(A) {
    const g = r(A).value;
    return x(g);
  }, B.wbg.__wbg_iterator_97f0c81209c6c35a = function() {
    return x(Symbol.iterator);
  }, B.wbg.__wbg_get_97b561fb56f034b5 = function() {
    return kA(function(A, g) {
      const C = Reflect.get(r(A), r(g));
      return x(C);
    }, arguments);
  }, B.wbg.__wbindgen_is_function = function(A) {
    return typeof r(A) == "function";
  }, B.wbg.__wbg_next_526fc47e980da008 = function(A) {
    const g = r(A).next;
    return x(g);
  }, B.wbg.__wbg_call_cb65541d95d71282 = function() {
    return kA(function(A, g) {
      const C = r(A).call(r(g));
      return x(C);
    }, arguments);
  }, B.wbg.__wbg_self_1ff1d729e9aae938 = function() {
    return kA(function() {
      const A = self.self;
      return x(A);
    }, arguments);
  }, B.wbg.__wbg_window_5f4faef6c12b79ec = function() {
    return kA(function() {
      const A = window.window;
      return x(A);
    }, arguments);
  }, B.wbg.__wbg_globalThis_1d39714405582d3c = function() {
    return kA(function() {
      const A = globalThis.globalThis;
      return x(A);
    }, arguments);
  }, B.wbg.__wbg_global_651f05c6a0944d1c = function() {
    return kA(function() {
      const A = global.global;
      return x(A);
    }, arguments);
  }, B.wbg.__wbg_isArray_4c24b343cb13cfb1 = function(A) {
    return Array.isArray(r(A));
  }, B.wbg.__wbg_instanceof_ArrayBuffer_39ac22089b74fddb = function(A) {
    let g;
    try {
      g = r(A) instanceof ArrayBuffer;
    } catch {
      g = !1;
    }
    return g;
  }, B.wbg.__wbg_call_01734de55d61e11d = function() {
    return kA(function(A, g, C) {
      const I = r(A).call(r(g), r(C));
      return x(I);
    }, arguments);
  }, B.wbg.__wbg_newwithbyteoffsetandlength_735ed5ea2ae07fe9 = function(A, g, C) {
    const I = new Int16Array(r(A), g >>> 0, C >>> 0);
    return x(I);
  }, B.wbg.__wbg_newwithbyteoffsetandlength_9f43b22ab631d1d6 = function(A, g, C) {
    const I = new Int32Array(r(A), g >>> 0, C >>> 0);
    return x(I);
  }, B.wbg.__wbg_newwithbyteoffsetandlength_6da8e527659b86aa = function(A, g, C) {
    const I = new Uint8Array(r(A), g >>> 0, C >>> 0);
    return x(I);
  }, B.wbg.__wbg_set_5cf90238115182c3 = function(A, g, C) {
    r(A).set(r(g), C >>> 0);
  }, B.wbg.__wbg_newwithbyteoffsetandlength_31ff1024ef0c63c7 = function(A, g, C) {
    const I = new Uint16Array(r(A), g >>> 0, C >>> 0);
    return x(I);
  }, B.wbg.__wbg_newwithbyteoffsetandlength_6df0e8c3efd2a5d3 = function(A, g, C) {
    const I = new Uint32Array(r(A), g >>> 0, C >>> 0);
    return x(I);
  }, B.wbg.__wbg_newwithbyteoffsetandlength_69193e31c844b792 = function(A, g, C) {
    const I = new Float32Array(r(A), g >>> 0, C >>> 0);
    return x(I);
  }, B.wbg.__wbg_instanceof_Uint8Array_d8d9cb2b8e8ac1d4 = function(A) {
    let g;
    try {
      g = r(A) instanceof Uint8Array;
    } catch {
      g = !1;
    }
    return g;
  }, B.wbg.__wbg_set_092e06b0f9d71865 = function() {
    return kA(function(A, g, C) {
      return Reflect.set(r(A), r(g), r(C));
    }, arguments);
  }, B.wbg.__wbindgen_jsval_loose_eq = function(A, g) {
    return r(A) == r(g);
  }, B.wbg.__wbindgen_bigint_from_i64 = function(A) {
    return x(A);
  }, B.wbg.__wbindgen_bigint_from_u64 = function(A) {
    const g = BigInt.asUintN(64, A);
    return x(g);
  }, B.wbg.__wbg_fromCodePoint_5aced4409dbb8969 = function() {
    return kA(function(A) {
      const g = String.fromCodePoint(A >>> 0);
      return x(g);
    }, arguments);
  }, B.wbg.__wbg_new_56693dbed0c32988 = function() {
    return x(/* @__PURE__ */ new Map());
  }, B.wbg.__wbg_isSafeInteger_bb8e18dd21c97288 = function(A) {
    return Number.isSafeInteger(r(A));
  }, B.wbg.__wbindgen_debug_string = function(A, g) {
    const C = sB(r(g)), I = JA(C, h.__wbindgen_malloc, h.__wbindgen_realloc), E = MA;
    f()[A / 4 + 1] = E, f()[A / 4 + 0] = I;
  }, B.wbg.__wbindgen_throw = function(A, g) {
    throw new Error(dA(A, g));
  }, B.wbg.__wbindgen_rethrow = function(A) {
    throw AA(A);
  }, B.wbg.__wbg_then_b2267541e2a73865 = function(A, g, C) {
    const I = r(A).then(r(g), r(C));
    return x(I);
  }, B.wbg.__wbg_resolve_53698b95aaf7fcf8 = function(A) {
    const g = Promise.resolve(r(A));
    return x(g);
  }, B.wbg.__wbg_then_f7e06ee3c11698eb = function(A, g) {
    const C = r(A).then(r(g));
    return x(C);
  }, B.wbg.__wbindgen_cb_drop = function(A) {
    const g = AA(A).original;
    return g.cnt-- == 1 ? (g.a = 0, !0) : !1;
  }, B.wbg.__wbg_byobRequest_08c18cee35def1f4 = function(A) {
    const g = r(A).byobRequest;
    return FA(g) ? 0 : x(g);
  }, B.wbg.__wbg_releaseLock_9ae075576f54bf0b = function() {
    return kA(function(A) {
      r(A).releaseLock();
    }, arguments);
  }, B.wbg.__wbg_close_e9110ca16e2567db = function(A) {
    r(A).close();
  }, B.wbg.__wbg_enqueue_d71a1a518e21f5c3 = function(A, g) {
    r(A).enqueue(r(g));
  }, B.wbg.__wbg_view_231340b0dd8a2484 = function(A) {
    const g = r(A).view;
    return FA(g) ? 0 : x(g);
  }, B.wbg.__wbg_byteLength_5299848ed3264181 = function(A) {
    return r(A).byteLength;
  }, B.wbg.__wbg_close_da7e6fb9d9851e5a = function(A) {
    r(A).close();
  }, B.wbg.__wbg_buffer_4e79326814bdd393 = function(A) {
    const g = r(A).buffer;
    return x(g);
  }, B.wbg.__wbg_byteOffset_b69b0a07afccce19 = function(A) {
    return r(A).byteOffset;
  }, B.wbg.__wbg_respond_8fadc5f5c9d95422 = function(A, g) {
    r(A).respond(g >>> 0);
  }, B.wbg.__wbg_canvas_2584a9e73f3f328d = function(A) {
    const g = r(A).canvas;
    return FA(g) ? 0 : x(g);
  }, B.wbg.__wbg_bindVertexArray_8863a216d7b0a339 = function(A, g) {
    r(A).bindVertexArray(r(g));
  }, B.wbg.__wbg_bufferData_21334671c4ba6004 = function(A, g, C, I) {
    r(A).bufferData(g >>> 0, r(C), I >>> 0);
  }, B.wbg.__wbg_bufferSubData_c472b93c9e272eac = function(A, g, C, I) {
    r(A).bufferSubData(g >>> 0, C, r(I));
  }, B.wbg.__wbg_readPixels_99fda83f6ca7ec72 = function() {
    return kA(function(A, g, C, I, E, o, t, s) {
      r(A).readPixels(g, C, I, E, o >>> 0, t >>> 0, r(s));
    }, arguments);
  }, B.wbg.__wbg_texImage2D_699c5d8e0d9ea28a = function() {
    return kA(function(A, g, C, I, E, o, t, s, M, N, k) {
      r(A).texImage2D(g >>> 0, C, I, E, o, t, s >>> 0, M >>> 0, N === 0 ? void 0 : VE(N, k));
    }, arguments);
  }, B.wbg.__wbg_bindBuffer_24f6010e273fa400 = function(A, g, C) {
    r(A).bindBuffer(g >>> 0, r(C));
  }, B.wbg.__wbg_bindFramebuffer_a9573e340dab20fe = function(A, g, C) {
    r(A).bindFramebuffer(g >>> 0, r(C));
  }, B.wbg.__wbg_bindTexture_92d6d7f8bff9531e = function(A, g, C) {
    r(A).bindTexture(g >>> 0, r(C));
  }, B.wbg.__wbg_blendEquation_12146cb96dc1bcd9 = function(A, g) {
    r(A).blendEquation(g >>> 0);
  }, B.wbg.__wbg_blendFuncSeparate_fbf93dee3e5ce456 = function(A, g, C, I, E) {
    r(A).blendFuncSeparate(g >>> 0, C >>> 0, I >>> 0, E >>> 0);
  }, B.wbg.__wbg_createBuffer_323425af422748ac = function(A) {
    const g = r(A).createBuffer();
    return FA(g) ? 0 : x(g);
  }, B.wbg.__wbg_createFramebuffer_1684a99697ac9563 = function(A) {
    const g = r(A).createFramebuffer();
    return FA(g) ? 0 : x(g);
  }, B.wbg.__wbg_createTexture_1bf4d6fec570124b = function(A) {
    const g = r(A).createTexture();
    return FA(g) ? 0 : x(g);
  }, B.wbg.__wbg_cullFace_6daa9f2aa42b4620 = function(A, g) {
    r(A).cullFace(g >>> 0);
  }, B.wbg.__wbg_deleteBuffer_2c09d03fa4b0bd08 = function(A, g) {
    r(A).deleteBuffer(r(g));
  }, B.wbg.__wbg_deleteFramebuffer_edd16bb8df6a8e0d = function(A, g) {
    r(A).deleteFramebuffer(r(g));
  }, B.wbg.__wbg_disable_e02106ca6c7002d6 = function(A, g) {
    r(A).disable(g >>> 0);
  }, B.wbg.__wbg_disableVertexAttribArray_6d57776c8f642f44 = function(A, g) {
    r(A).disableVertexAttribArray(g >>> 0);
  }, B.wbg.__wbg_enable_195891416c520019 = function(A, g) {
    r(A).enable(g >>> 0);
  }, B.wbg.__wbg_framebufferTexture2D_e88fcbd7f8523bb8 = function(A, g, C, I, E, o) {
    r(A).framebufferTexture2D(g >>> 0, C >>> 0, I >>> 0, r(E), o);
  }, B.wbg.__wbg_getProgramParameter_35522a0bfdfaad27 = function(A, g, C) {
    const I = r(A).getProgramParameter(r(g), C >>> 0);
    return x(I);
  }, B.wbg.__wbg_texParameteri_85dad939f62a15aa = function(A, g, C, I) {
    r(A).texParameteri(g >>> 0, C >>> 0, I);
  }, B.wbg.__wbg_uniform1f_88379f4e2630bc66 = function(A, g, C) {
    r(A).uniform1f(r(g), C);
  }, B.wbg.__wbg_uniform1i_d2e61a6a43889648 = function(A, g, C) {
    r(A).uniform1i(r(g), C);
  }, B.wbg.__wbg_uniform2f_b6e484a1302ea599 = function(A, g, C, I) {
    r(A).uniform2f(r(g), C, I);
  }, B.wbg.__wbg_useProgram_3683cf6f60939dcd = function(A, g) {
    r(A).useProgram(r(g));
  }, B.wbg.__wbg_viewport_fad1ce9e18f741c0 = function(A, g, C, I, E) {
    r(A).viewport(g, C, I, E);
  }, B.wbg.__wbg_document_f7ace2b956f30a4f = function(A) {
    const g = r(A).document;
    return FA(g) ? 0 : x(g);
  }, B.wbg.__wbg_fetch_336b6f0cb426b46e = function(A, g) {
    const C = r(A).fetch(r(g));
    return x(C);
  }, B.wbg.__wbg_getwithindex_a41abfefef0f0eb7 = function(A, g) {
    const C = r(A)[g >>> 0];
    return FA(C) ? 0 : x(C);
  }, B.wbg.__wbg_newwithstrandinit_cad5cd6038c7ff5d = function() {
    return kA(function(A, g, C) {
      const I = new Request(dA(A, g), r(C));
      return x(I);
    }, arguments);
  }, B.wbg.__wbg_instanceof_HtmlCanvasElement_da5f9efa0688cf6d = function(A) {
    let g;
    try {
      g = r(A) instanceof HTMLCanvasElement;
    } catch {
      g = !1;
    }
    return g;
  }, B.wbg.__wbg_getContext_7c5944ea807bf5d3 = function() {
    return kA(function(A, g, C) {
      const I = r(A).getContext(dA(g, C));
      return FA(I) ? 0 : x(I);
    }, arguments);
  }, B.wbg.__wbg_instanceof_Window_9029196b662bc42a = function(A) {
    let g;
    try {
      g = r(A) instanceof Window;
    } catch {
      g = !1;
    }
    return g;
  }, B.wbg.__wbg_setcrossOrigin_99cf79991d1a3c2d = function(A, g, C) {
    r(A).crossOrigin = g === 0 ? void 0 : dA(g, C);
  }, B.wbg.__wbg_new_6f9cb260fad32a20 = function() {
    return kA(function() {
      const A = new Image();
      return x(A);
    }, arguments);
  }, B.wbg.__wbg_arrayBuffer_288fb3538806e85c = function() {
    return kA(function(A) {
      const g = r(A).arrayBuffer();
      return x(g);
    }, arguments);
  }, B.wbg.__wbg_setProperty_b95ef63ab852879e = function() {
    return kA(function(A, g, C, I, E) {
      r(A).setProperty(dA(g, C), dA(I, E));
    }, arguments);
  }, B.wbg.__wbindgen_closure_wrapper748 = function(A, g, C) {
    const I = oQ(A, g, 76, mE);
    return x(I);
  }, B.wbg.__wbindgen_closure_wrapper1367 = function(A, g, C) {
    const I = oQ(A, g, 146, OE);
    return x(I);
  }, B.wbg.__wbindgen_closure_wrapper1692 = function(A, g, C) {
    const I = vE(A, g, 76, bE);
    return x(I);
  }, B;
}
function UQ(B, A) {
  return h = B.exports, MB.__wbindgen_wasm_module = A, yI = null, GI = null, hI = null, MI = null, cI = null, h;
}
function Di(B) {
  if (h !== void 0)
    return h;
  const A = JQ();
  B instanceof WebAssembly.Module || (B = new WebAssembly.Module(B));
  const g = new WebAssembly.Instance(B, A);
  return UQ(g, B);
}
async function MB(B) {
  if (h !== void 0)
    return h;
  const A = JQ();
  (typeof B == "string" || typeof Request == "function" && B instanceof Request || typeof URL == "function" && B instanceof URL) && (B = fetch(B));
  const { instance: g, module: C } = await oi(await B, A);
  return UQ(g, C);
}
const ai = /* @__PURE__ */ Object.freeze(/* @__PURE__ */ Object.defineProperty({
  __proto__: null,
  AngleSerializeFmt: XE,
  BlendCfg: nI,
  BlendFactor: $E,
  BlendFunc: Ai,
  CenteredFoV: gi,
  ColorRGB: FI,
  ColorRGBA: vg,
  CooSystem: zE,
  GridCfg: Ii,
  ImageExt: _E,
  ImageMetadata: vI,
  IntoUnderlyingByteSource: Bi,
  IntoUnderlyingSink: Qi,
  IntoUnderlyingSource: Ci,
  MOC: ng,
  PipeOptions: Ei,
  QueuingStrategy: ii,
  ReadableStreamGetReaderOptions: bI,
  TransferFunction: PE,
  WebClient: mI,
  default: MB,
  initSync: Di
}, Symbol.toStringTag, { value: "Module" }));
let SA = {};
SA.aladin = function(B, A) {
  return new VA(O(B)[0], A);
};
SA.source = function(B, A, g, C) {
  return new GB(B, A, g, C);
};
SA.marker = function(B, A, g, C) {
  return g = g || {}, g.marker = !0, SA.source(B, A, C, g);
};
SA.polygon = function(B, A) {
  const g = B.length;
  if (g < 3)
    throw "Cannot define a polygon from less than 3 vertices";
  const C = g - 1;
  return B[0][0] == B[C][0] && B[0][1] == B[C][1] && B.pop(), A = A || {}, A.closed = !0, new eB(B, A);
};
SA.polyline = function(B, A) {
  return new eB(B, A);
};
SA.circle = function(B, A, g, C) {
  return new hQ([B, A], g, C);
};
SA.ellipse = function(B, A, g, C, I, E) {
  return new MQ([B, A], g, C, I, E);
};
SA.graphicOverlay = function(B) {
  return new gI(B);
};
SA.catalog = function(B) {
  return new bg(B);
};
SA.catalogHiPS = function(B, A) {
  return new NQ(B, null, null, A);
};
SA.footprintsFromSTCS = function(B, A) {
  var g = gI.parseSTCS(B, A);
  return g;
};
SA.MOCFromURL = function(B, A, g) {
  var C = new RQ(A);
  return C.dataFromFITSURL(B, g), C;
};
SA.MOCFromJSON = function(B, A) {
  var g = new RQ(A);
  return g.dataFromJSON(B), g;
};
SA.catalogFromURL = function(B, A, g, C, I) {
  var E = SA.catalog(A);
  const o = function(t, s, M) {
    E.setFields(M), E.isObsCore() && (E.name = "ObsCore:" + B), rQ.handleActions(E), E.addFootprints(s), E.addSources(t), g && g(t);
  };
  return I !== void 0 ? bg.parseVOTable(
    B,
    o,
    C,
    E.maxNbSources,
    I,
    E.raField,
    E.decField
  ) : bg.parseVOTable(
    B,
    o,
    () => {
      console.log("error cors"), bg.parseVOTable(
        B,
        o,
        C,
        E.maxNbSources,
        !0,
        E.raField,
        E.decField
      );
    },
    E.maxNbSources,
    !1,
    E.raField,
    E.decField
  ), E;
};
SA.catalogFromSimbad = function(B, A, g, C, I) {
  g = g || {}, "name" in g || (g.name = "Simbad");
  var E = zg.buildSimbadCSURL(B, A);
  return SA.catalogFromURL(E, g, C, I, !1);
};
SA.catalogFromNED = function(B, A, g, C, I) {
  g = g || {}, "name" in g || (g.name = "NED");
  var E;
  if (B && typeof B == "object")
    "ra" in B && "dec" in B && (E = zg.buildNEDPositionCSURL(B.ra, B.dec, A));
  else {
    var o = /[a-zA-Z]/.test(B);
    if (o)
      E = zg.buildNEDObjectCSURL(B, A);
    else {
      var t = new $A();
      t.parse(B), E = zg.buildNEDPositionCSURL(t.lon, t.lat, A);
    }
  }
  return SA.catalogFromURL(E, g, C, I, !0);
};
SA.catalogFromVizieR = function(B, A, g, C, I, E) {
  C = C || {}, "name" in C || (C.name = "VizieR:" + B);
  var o = zg.buildVizieRCSURL(B, A, g, C);
  return SA.catalogFromURL(o, C, I, E, !1);
};
SA.catalogFromSkyBot = function(B, A, g, C, I, E, o, t) {
  I = I || {}, E = E || {}, "name" in E || (E.name = "SkyBot");
  var s = zg.buildSkyBotCSURL(B, A, g, C, I);
  return SA.catalogFromURL(s, E, o, t, !1);
};
SA.hipsDefinitionFromURL = function(B, A) {
  kQ.fromURL(B, A);
};
SA.getAvailableListOfColormaps = function() {
  return OI.COLORMAPS;
};
SA.init = (async () => {
  const B = document.createElement("canvas").getContext("webgl2");
  if (await MB(), B)
    VA.wasmLibs.core = ai;
  else
    throw "WebGL2 not supported by your browser";
})();
const zA = SA;
export {
  zA as default
};