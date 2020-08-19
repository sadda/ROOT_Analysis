random_points(n::Int, pars::Pars) = pars.x_min .+ (pars.x_max .- pars.x_min).*rand(pars.d, n)


get_file_name(folder_name, δ, m) = joinpath(folder_name, "Results_delta=" * string(δ) * "_m=" * string(m) * ".bson")


create_directory(folder_name) = !isdir(folder_name) && mkdir(folder_name)


using Printf


print_formatted(fmt, args...) = @eval @printf($fmt, $(args...))


function table_to_tex(data, format; header_l=[], header_t=[], caption="", label="", alignment=[])

    if ndims(data) == 1
        data_aux = repeat(data, 1, 1)
    else
        data_aux = data
    end

    if ndims(data_aux) != 2
        error("The dimension of data must be 1 or 2")
    end

    n_row, n_col = size(data_aux)

    if typeof(format) <: String
        format_all = [format for i in 1:n_col]
    elseif size(format) == 1
        format_all = [format[1] for i in 1:n_col]
    end

    #Beginning
    @printf("\n\n")
    @printf("\\begin{table}[!ht]\n")
    @printf("\\caption{%s}\n", caption)
    @printf("\\label{%s}\n", label)
    @printf("\\centering\n")
    if isempty(alignment)
        if isempty(header_l)
            alignment = "@{}" * repeat("l", n_col) * "@{}"
        else
            alignment = "@{}" * repeat("l", n_col+1) * "@{}"
        end
    end
    @printf("\\begin{tabular}{%s}\n", alignment)
    @printf("\\toprule\n")

    # Header
    if !isempty(header_t)
        if typeof(header_t) <: String
            @printf("%s\n", header_t)
        else
            if !isempty(header_l)
                @printf(" & ")
            end
            for i in 1:length(header_t)
                if typeof(header_t[i]) <: Number
                    @printf("\$%f\$", header_t[i])
                else
                    @printf("%s", header_t[i])
                end
                if i < n_col
                    @printf(" & ")
                end
            end
        end
        @printf(" \\\\\n\\midrule\n")
    end

    # Content
    for i=1:n_row
        if !isempty(header_l)
            try
                @printf("%s", header_l[i])
            catch
                @printf("\$%f\$", header_l[i])
            end
            @printf(" & ")
        end

        for j in 1:n_col
            value       = data[i,j]
            format_type = format_all[j]
            if format_type[end] == 'e'

            elseif format_type[end] == 'p'
                format_type = format_type[1:end-1] * 'f'
                @printf("\$")
                print_formatted(format_type, 100*value)
                @printf("\\%s\$", "%")
            elseif format_type[end] == 's'

            else
                @printf("\$")
                print_formatted(format_type, value)
                @printf("\$")
            end
            if j < n_col
                @printf(" & ")
            else
                @printf(" \\\\\n")
            end
        end
    end

    # End
    @printf("\\bottomrule\n")
    @printf("\\end{tabular}\n")
    @printf("\\end{table}\n")
    @printf("\n\n")
end


alignment_lr(x::AbstractArray) = "@{}l" * repeat('r', length(x)) * "@{}"
