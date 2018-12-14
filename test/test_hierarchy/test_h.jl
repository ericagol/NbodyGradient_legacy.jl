# Sample of how to run the arbirary hierarchy generator "setup_hierarchy.jl"
# This example uses the file "test.txt", which is a sample initial condition
# file.

#The function also takes an array if you wish to use it in the terminal.
# They array should appear as: name = [x, "a,b,c,..."], the first element of the
# array should be an Int64 and the second a String of integers separated by
#commas. x it the number of bodies, and a,b,c,... are the number of binaries on
# a level of the hierarchy (See "test.txt" for reference).


include("setup_hierarchy.jl")

# Example of using the file method and outputs to STDOUT
file = "test.txt"
h = hierarchy(file)
print("From IC file: ""test.txt"" \n")
Base.showarray(STDOUT,h,false)
print("\n")

# Example of same hierarchy, using the direct array input method
file = [7, "1,2,2,1"]
h = hierarchy(file)
print("From direct array input ([7, ""1,2,2,1""]) \n")
Base.showarray(STDOUT,h,false)
print("\n")

# Example of hierarchy 1 in "hierarchies.txt"
h = hierarchy([6,"2,2,1"])
print("From direct array input ([6, ""2,2,1""]) \n")
Base.showarray(STDOUT,h,false)
print("\n")

# Example of hierarchy 2
h = hierarchy([5,"1,1,1,1"])
print("From direct array input ([5, ""1,1,1,1""]) \n")
Base.showarray(STDOUT,h,false)
print("\n")

#Example of hierarchy 3
h = hierarchy([8,"4,2,1"])
print("From direct array input ([6, ""2,2,1""]) \n")
Base.showarray(STDOUT,h,false)
print("\n")
