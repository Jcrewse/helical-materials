# Path to VMD application
VMD_PATH="/Applications/VMD\ 1.9.4a57-x86_64-Rev12.app"

# TCL script to run from VMD
# TCL_SCRIPT=$(cat << 'EOF'
# # Set the path to the folder containing your .cube files
# set folder_path "/Users/jcrewse/Research/Helical"

# # Set the file prefix for your .cube files
# set file_prefix "H2-Nphi-1"

# # Set the range of unique identifiers for your files
# set start_id 0
# set end_id 1

# # Clear previously loaded molecules
# mol delete all

# # Iterate through the files
# for {set i $start_id} {$i <= $end_id} {incr i} {
#     # Load the .cube file
#     set file_name "${folder_path}/${file_prefix}_${i}_RePsi.cube"
#     puts $file_name
#     mol new $file_name type cube first 0 last -1 step 1 waitfor all

#     # Change the graphical representation to CPK
#     mol delrep 0 top
#     mol representation CPK
#     mol addrep top

#     # Add the isosurface representation
#     mol representation Isosurface
#     mol addrep top
#     set rep_id [expr {[molinfo top get numreps] - 1}]
#     mol modstyle $rep_id top "Isosurface 0.1 0 0 0 0"
#     mol modmaterial $rep_id top Transparent
# }

# # Center and zoom to fit all molecules
# mol modselect 0 top "all"
# mol modstyle 0 top "CPK"
# mol modcolor 0 top "ColorID 6"

# # Rotate and zoom out
# rotate y by 40
# scale by 0.56
# )

echo "$TCL_SCRIPT" | $VMD_PATH -dispdev text -args
echo $VMD_PATH -dispdev text -args