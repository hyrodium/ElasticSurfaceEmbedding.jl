function SlackString(str::String)
    str2 = "{\"text\":\""*str*"\"}"
    cmd = `curl -X POST -H 'Content-type: application/json' --data $str2 $IWhU`
    # println(cmd)
    run(cmd; wait=false)
    return nothing
end

function SlackDict(dic::Dict)
    DIC = ""
    for key in sort(collect(keys(dic)))
        DIC = DIC*key*" : "*string(dic[key])*"\n"
    end
    DIC = "{\"text\":\""*DIC*"\"}"
    cmd = `curl -X POST -H 'Content-type: application/json' --data $DIC $IWhU`
    # println(cmd)
    run(cmd; wait=false)
    return nothing
end

function SlackFile(filename; comment="")
    cmd = `curl -F file=@$filename -F "initial_comment=$comment" -F channels=$ChID -H "Authorization: Bearer $OAAT" https://slack.com/api/files.upload`
    # println(cmd)
    run(cmd; wait=false)
    return nothing
end
