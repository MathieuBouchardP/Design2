app = AppClass;
app.state = 0;

function simple_a(var_com)
    for i = 1:10
        disp(var_com.state)
        pause(2)
    end
end

simple_a(app);
