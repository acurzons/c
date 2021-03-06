ds9 ../rxj1713_pearson_suz_mHI.fits \
-cmap bb \
-colorbar yes \
-colorbar vertical \
-colorbar fontsize 16 \
-colorbar size 15 \
-colorbar space value \
-colorbar ticks 13 \
-contour load ../../DS9/hess_7_16_27.ctr green 2 no \
-crop +17:13:44.0514 -39:47:13.886 98 93 wcs fk5 arcmin \
-scale limits -1 1 \
-grid yes \
-grid system wcs \
-grid sky fk5 \
-grid skyformat degrees \
-grid type publication \
-grid type axes exterior \
-grid type numerics interior \
-grid view grid no \
-grid view axes yes \
-grid view title no \
-grid format1 d.1 \
-grid format2 d.1 \
-grid numerics color white \
-grid numerics gap2 6 \
-grid grid gap1 0.6 \
-grid grid gap2 0.6 \
-grid labels gap2 3 \
-zoom to fit \
-width 428 \
-height 414 \
-pan to 17:14:04 -39:47:43 wcs fk5 \
-saveimage eps rxj1713_pearson_suz_mHI.ps

