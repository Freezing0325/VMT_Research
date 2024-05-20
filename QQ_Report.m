function QQ_Report(qqnum,message)
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