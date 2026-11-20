import numpy as np
from PIL import Image
def calib(fn, vcol, hrow):
    im=np.asarray(Image.open(f'figs/{fn}_flipped.png').convert('L')).astype(int)
    # vertical-axis major ticks: rows where 4 px right of axis are dark (=102) and next are white
    vt=[r for r in range(im.shape[0]) if (im[r,vcol+1:vcol+5]<=110).all() and im[r,vcol+6]>=250]
    # horizontal-axis major ticks: cols where 4 px above axis are dark
    ht=[c for c in range(im.shape[1]) if (im[hrow-5:hrow-1,c]<=110).all() and im[hrow-7,c]>=250]
    return vt, ht
print('lambda', calib('lambda', 23, 609))
print('lambda2', calib('lambda2', 40, 512))
