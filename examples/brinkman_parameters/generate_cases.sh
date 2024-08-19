# omg Tim you're going to read this eventually and this it's such a stupid way of doing this hahahaha

# Mesh		4
# dt			10
# Re			25
# Chi			69
# implcit	72
# radius		76




# Study
declare -a mesh_list=("2" "3" "4")
declare -a dt_list=("2.5e-3" "1.6e-3" "1.25e-3")
declare -a Re_list=("200")
declare -a chi_list=("1" "10" "100" "1000" "1000000")
declare -a implicit_list=("true")
declare -a radius_list=("0.05")

big_name="YOFAM"
mkdir $big_name

for mesh in "${mesh_list[@]}"
do
for Re in "${Re_list[@]}"
do
for chi in "${chi_list[@]}"
do
for implicit in "${implicit_list[@]}"
do
for radius in "${radius_list[@]}"
do
  name=$big_name/mesh${mesh}_Re${Re}_chi${chi}_radius${radius}_$implicit
  cp -r default_case $name
	# now all the replacements

	new_line='4s#.*#"mesh_file":"../../MESHES/immersed_M'$mesh'.nmsh",#'
	sed -i $new_line $name/cylinder.case
	new_line='25s#.*#"Re":'$Re',#'
	sed -i $new_line $name/cylinder.case
	new_line='69s#.*#'$chi',#'
	sed -i $new_line $name/cylinder.case
	new_line='72s#.*#"implicit":'$implicit',#'
	sed -i $new_line $name/cylinder.case
	new_line='76s#.*#"radius":'$radius',#'
	sed -i $new_line $name/cylinder.case

	# time step is a bit strange
	case $mesh in
	"2")
	dt="2.5e-3"
	;;
	"3")
	dt="1.6e-3"
	;;
	"4")
	dt="1.25e-3"
	;;
	esac
	new_line='10s#.*#"timestep":'$dt',#'
	sed -i $new_line $name/cylinder.case

done
done
done
done
done


