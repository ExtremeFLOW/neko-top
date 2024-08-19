# omg Tim you're going to read this eventually and this it's such a stupid way of doing this hahahaha

# Mesh		4
# dt			10
# Re			25
# Chi			69
# implcit	72
# radius		76




# Study
declare -a mesh_list=("2" "3" "4")
declare -a Re_list=("200")

big_name="MESHED"
mkdir $big_name

for mesh in "${mesh_list[@]}"
do
for Re in "${Re_list[@]}"
do
  name=$big_name/mesh${mesh}_Re${Re}
  cp -r default_case_meshed $name
	# now all the replacements

	new_line='4s#.*#"mesh_file":"../../../../data/brinkman_parameters/meshed_M'$mesh'.nmsh",#'
	sed -i $new_line $name/cylinder.case
	new_line='25s#.*#"Re":'$Re',#'
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


