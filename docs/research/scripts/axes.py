import numpy as np
from PIL import Image
for fn in ['lambda','lambda2']:
    im = np.asarray(Image.open(f'figs/{fn}_flipped.png').convert('RGB')).astype(int)
    H,W,_=im.shape
    grey = (np.abs(im[:,:,0]-im[:,:,1])<12)&(np.abs(im[:,:,1]-im[:,:,2])<12)&(im[:,:,0]<200)&(im[:,:,0]>60)
    colsum = grey.sum(axis=0); rowsum=grey.sum(axis=1)
    vcol = int(np.argmax(colsum[:300])); hrow=int(np.argmax(rowsum))
    print(fn,'vertical axis col',vcol,'count',colsum[vcol],' horizontal axis row',hrow,'count',rowsum[hrow])
    # ticks on vertical axis: grey pixels in columns vcol+2..vcol+8
    tick = grey[:, vcol+2:vcol+9].sum(axis=1)
    rows = np.where(tick>=5)[0]
    # group
    groups=[]; 
    for r in rows:
        if groups and r-groups[-1][-1]<=1: groups[-1].append(r)
        else: groups.append([r])
    print('  long ticks (rows):',[round(np.mean(g),1) for g in groups])
    tick2 = grey[hrow-8:hrow-1,:].sum(axis=0)
    cols=np.where(tick2>=5)[0]; groups=[]
    for c in cols:
        if groups and c-groups[-1][-1]<=1: groups[-1].append(c)
        else: groups.append([c])
    print('  long ticks on h-axis (cols):',[round(np.mean(g),1) for g in groups])
