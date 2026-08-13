# Addendum (2026-08-13): what the unlabeled 55 percent is

A follow-up measurement (results/voidct1218) checked the CT under the
m7-positive voxels, using the masked canonical volume with the same
timestamp as the prediction. 46.5 percent of ALL m7-positive voxels sit in
CT chunks that do not exist in the masked volume (85.4 percent inside
void-dominated windows, exactly 0.0 in fully-labeled controls), and the m7
fill value is 0, so these positives are real stored predictions over
regions the masked reconstruction leaves empty. The coverage number in
this directory therefore splits: most of the unlabeled predicted surface
is prediction beyond masked-CT support, and the genuine label gap is
small and concentrated in the lowest z band, where void windows show real
sheet-intensity CT. Consumers thresholding the m7 surface should
intersect it with CT support first.
