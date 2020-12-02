# --
# Given a T value, and single component data, compute the psi value for each component in the mixture.
# --
function compute_psi(T::Float64, data_array::Array{Float64,2})::Array{Float64,1}

    # ok, so we need to compute the psi_* value -
    R = 8.314 # J/mol-K

    # what is the size of the data array?
    # components are on the rows
    # DH = col 1
    # Tm = col 2
    (number_of_components, number_of_cols) = size(data_array)

    # initialize -
    psi_vector = zeros(number_of_components)
    
    # main -
    for component_index = 1:number_of_components
        
        # get data for component_index =
        𝝙H = data_array[component_index,1]
        Tm = data_array[component_index,2]

        # compute psi for component_index -
        psi_value = -(𝝙H/R)*(1/T - 1/Tm)
    
        # package -
        psi_vector[component_index] = exp(psi_value)
    end

    # return -
    return psi_vector
end

# --
# Given an array of T values, and single component data, compute the psi value for each component in the mixture.
# Calls helper function for each value of T
# --
function compute_psi(T_vector::Array{Float64,1},data_array::Array{Float64,2})::Array{Float64,2}

    # initialize -
    number_of_temperature_points = length(T_vector)
    psi_array = Array{Float64,2}(undef,number_of_temperature_points,2)

    # main -
    for (index,T_value) in enumerate(T_vector)
        
        # compute psi_1 and psi_2 -
        psi_vec_tmp = compute_psi(T_value, data_array)
        
        # package -
        psi_array[index,1] = psi_vec_tmp[1]
        psi_array[index,2] = psi_vec_tmp[2]
    end

    # return -
    return psi_array
end

# --
# Computes the composition x1 and z1 given psi value array. Uses SVN eqn 15.20 and 15.21 (derived from SLE matching, ideal phases)
# --
function compute_composition(psi_array::Array{Float64,2})::Array{Float64,2}

    # initialize -
    (number_of_rows,number_of_cols) = size(psi_array)
    composition_array = Array{Float64,2}(undef,number_of_rows, number_of_cols)
    
    # main -
    # rows = temperatures
    # cols = species 
    for row_index = 1:number_of_rows
        
        # get the 𝜓 values -
        psi_1_value = psi_array[row_index,1]
        psi_2_value = psi_array[row_index,2]

        # compute the compostion -
        # use SVN eqn 15.20 and 15.21 (derived from SLE matching, ideal phases)
        x1 = (psi_1_value)*(1-psi_2_value)/(psi_1_value - psi_2_value)
        z1 = (1-psi_2_value)/(psi_1_value - psi_2_value)

        # package -
        composition_array[row_index,1] = x1
        composition_array[row_index,2] = z1
    end

    # return -
    return composition_array
end