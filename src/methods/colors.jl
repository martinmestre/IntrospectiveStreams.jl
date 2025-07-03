function wongcolors_ext()
    return [
        RGB(0/255, 114/255, 178/255), # blue
        RGB(230/255, 159/255, 0/255), # orange
        RGB(0/255, 158/255, 115/255), # green
        RGB(204/255, 121/255, 167/255), # reddish purple (Slightly desaturated pink)
        RGB(86/255, 180/255, 233/255), # sky blue
        RGB(213/255, 94/255, 0/255), # vermillion
        RGB(240/255, 228/255, 66/255), # yellow
        RGB(227/255, 66/255, 52/255), # real vermillion
        RGB(148/255, 103/255, 189/255), # Lavander purple (Slightly desaturated violet)
        RGB(188/255, 189/255, 34/255), # Olive Yellow  (Strong yellow)
        RGB(128/255, 128/255, 0), #dark yellow (olive tone)
        RGB(190/155, 0, 77/155) # Cherry Red
    ]
end
function paulcolors()
    # https://www.nceas.ucsb.edu/sites/default/files/2022-06/Colorblind%20Safe%20Color%20Schemes.pdf
    return [
        RGB(221/255, 221/255, 221/255), # grey
        RGB(46/255, 37/255, 133/255), # blue
        RGB(51/255, 117/255, 56/255), #  green
        RGB(93/255, 168/255, 153/255), # casheza :)
        RGB(148/255, 203/255, 236/255), # celeste
        RGB(220/255, 205/255, 125/255), # amarillo
        RGB(194/255,106/255, 119/255), # rosa
        RGB(159/255, 74/255, 150/255), # violeta
        RGB(126/255, 41/255, 84/255) # morado
    ]
end
function wongcolors()
    return [
        RGB(0/255, 114/255, 178/255), # blue
        RGB(230/255, 159/255, 0/255), # orange
        RGB(0/255, 158/255, 115/255), # green
        RGB(204/255, 121/255, 167/255), # reddish purple
        RGB(86/255, 180/255, 233/255), # sky blue
        RGB(213/255, 94/255, 0/255), # vermillion
        RGB(240/255, 228/255, 66/255), # yellow
    ]
end

export wongcolors, wongcolors_ext, paulcolors