import fitz, glob, os
for f in sorted(glob.glob('src/*.pdf')):
    doc = fitz.open(f)
    pg = doc[0]
    dr = pg.get_drawings()
    imgs = pg.get_images()
    txt = pg.get_text()[:200].replace('\n',' | ')
    npts = sum(len(d['items']) for d in dr)
    print(f"{os.path.basename(f):18s} drawings={len(dr):5d} items={npts:7d} images={len(imgs)} rect={pg.rect} text='{txt}'")
