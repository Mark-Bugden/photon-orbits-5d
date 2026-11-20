import pymupdf, glob, os
for f in sorted(glob.glob('src/*.pdf')):
    doc = pymupdf.open(f)
    for img in doc[0].get_images():
        xref = img[0]
        pix = pymupdf.Pixmap(doc, xref)
        if pix.n - pix.alpha >= 4: pix = pymupdf.Pixmap(pymupdf.csRGB, pix)
        out = 'figs/' + os.path.basename(f).replace('.pdf', '.png')
        pix.save(out)
        print(out, pix.width, pix.height)
