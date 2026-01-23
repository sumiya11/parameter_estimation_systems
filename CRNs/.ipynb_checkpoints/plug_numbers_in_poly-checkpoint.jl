using Nemo

function make_dummy_system(sys_y, sys_x, main_vars, main_vars_val, diff_vars)
    sys_current = copy(sys_x)
    vars_current = copy(main_vars)
    vals_current = copy(main_vars_val)
    for i in 1:10^4
      sys_current = map(eq -> evaluate(eq, vars_current, vals_current), sys_current)
      sys_current = filter(!iszero, sys_current)
      for j in 1:length(sys_current)
        if leading_term(sys_current[j]) in setdiff(diff_vars, vars_current)
          if !(length(sys_current[j]) in (1,2)) && continue end
          if !issubset(vars(sys_current[j]), vcat(leading_term(sys_current[j]), vars_current)) && continue end
          if length(sys_current[j]) == 1
            new_var = leading_term(sys_current[j])
            new_val = zero(trailing_coefficient(sys_current[j]))
          else
            new_var = leading_term(sys_current[j])
            new_val = -trailing_coefficient(sys_current[j])
          end
          push!(vars_current, new_var)
          push!(vals_current, new_val)
        end
      end
      if Set(vars_current) == Set(vcat(main_vars, diff_vars))
        break
      end
    end

    @assert isempty(setdiff(diff_vars, vars_current))
    sys_y_dummy = [
      (f = eq - trailing_coefficient(eq); f + evaluate(-f, vars_current, vals_current))
      for eq in sys_y
    ]
    
    sys = vcat(sys_y_dummy, sys_x)
    @assert all(iszero, map(eq -> evaluate(eq, vars_current, vals_current), sys))
  
    sys, vars_current, vals_current
end

function plug_numbers_in_poly(sys, states, parameters)
  present_vars = vcat(states, parameters) # Base.union(map(vars, sys_x)...)
  present_vars = filter(x -> any(f -> degree(f, x) > 0, sys), present_vars)

  sys_y = filter(s -> !issubset(vars(s), present_vars), sys)
  sys_x = setdiff(sys, sys_y)
  
  main_vars = filter(x -> !endswith(string(x), r"_[1-9]"), present_vars)
  diff_vars = setdiff(present_vars, main_vars)

  # println("main_vars = \n", main_vars)
  # println("diff_vars = \n", diff_vars)
  # @info "" sys_y sys_x

  vals_dummy = [QQ(rand(1:100)) for i in 1:length(main_vars)]
  
  sys_dummy, vars_dummy, vals_dummy = make_dummy_system(sys_y, sys_x, main_vars, vals_dummy, diff_vars)
  
  sys_dummy, Dict(vars_dummy .=> vals_dummy)
end

begin
  R, (x_0, x_1, a, y_0, y_1) = QQ["x_0", "x_1", "a", "y_0", "y_1"];
  sys = [x_1 - x_0*a, y_0 - x_0^2, y_1 - 2*x_1*x_0]
  sys_num, subs = plug_numbers_in_poly(sys, [x_0, x_1], [a])

  R, (x_0, x_1, a, y_0, y_1) = QQ["x_0", "x_1", "a", "y_0", "y_1"];
  sys = [x_1 - 13*a + a^2, y_0 - 13^2 - a^2]
  sys_num, subs = plug_numbers_in_poly(sys, [x_0, x_1], [a])
end