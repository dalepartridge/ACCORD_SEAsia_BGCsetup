
for i in bdyfiles/accord*nc; do
  ncks -A bdy_gdept.nc $i
done

