function QQ_Report(qqnum,message)
%QQ_Report      利用QQ汇报一条消息
%   
%   qqnum       目标QQ号，即消息发送的对象
%   message     消息字符串


    if (nargin==1)
        message = '程序监控测试';
    end
    vb = actxserver('wscript.shell');
    eval(['!start tencent://Message/?Uin=',qqnum]);
    pause(1);
    clipboard('copy',message);
    pause(0.1);
    vb.SendKeys('^a');
    vb.SendKeys('{Delete}');
    vb.SendKeys('^v');
    pause(0.1);
    vb.SendKeys('^{ENTER}');
end