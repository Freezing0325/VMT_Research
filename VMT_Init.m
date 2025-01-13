function VMT_Init(CalMode)
    global a Normal_h CompPlus_h CompMinus_h OutputE Output_h Output_a;
        if (~exist('CalMode', 'var'))
            CalMode = 0;
        end
    switch CalMode
        case 0
        
            a = 35;
            Normal_h = 21.354;
            CompMinus_h = 11.368;
            CompPlus_h = 31.354;
            OutputE = 30.957;
            Output_h = 5;
            Output_a = 50;
            fprintf('实验模式参数设置完毕。\n');
        case 1
            a = 2;
            Normal_h = 1;
            CompMinus_h = 0.75;
            CompPlus_h = 1.25;
            OutputE = 600;
            Output_h = 0.125;
            Output_a = 2.5;
            fprintf('有限元模式参数设置完毕。\n');
    end

end