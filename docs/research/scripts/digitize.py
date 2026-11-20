import numpy as np
from PIL import Image
cols = {'x':(0x5E,0x81,0xB5),'theta':(0xE1,0x9C,0x24),'phi':(0x8F,0xB0,0x32),'psi':(0xEB,0x62,0x35)}
for fn in ['lambda','lambda2']:
    im = np.asarray(Image.open(f'figs/{fn}.png').convert('RGB')).astype(int)[::-1]  # flip vertically
    Image.fromarray(im.astype(np.uint8)).save(f'figs/{fn}_flipped.png')
    H,W,_ = im.shape
    print(fn, W, H)
    curves={}
    for k,c in cols.items():
        d = np.abs(im - np.array(c)).sum(axis=2)
        mask = d < 60
        xs=[];ys=[]
        for col in range(W-260):  # exclude legend at right
            rows=np.where(mask[:,col])[0]
            if len(rows): xs.append(col); ys.append(rows.mean())
        curves[k]=(np.array(xs),np.array(ys))
        print(f"  {k}: columns {xs[0] if xs else None}-{xs[-1] if xs else None}, row at start {ys[0] if ys else None:.1f}, row range {min(ys):.1f}-{max(ys):.1f}")
    np.save(f'figs/{fn}_curves.npy', curves, allow_pickle=True)
