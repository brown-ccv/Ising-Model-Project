# Function to write x-y data to a text file
function write_xy_data_to_file(filename, x_data, y_data)
    open(filename, "w") do file
        for (x, y) in zip(x_data, y_data)
            write(file, string(x) * "\t" * string(y) * "\n")
        end
    end
end

# Generate x-y data
x_data = collect(1:10)
y_data = x_data .^ 2

# Write x-y data to a text file
output_file = "xy_output.txt"
write_xy_data_to_file(output_file, x_data, y_data)

println("X-Y data has been written to $output_file")
