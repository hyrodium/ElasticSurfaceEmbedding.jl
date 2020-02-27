function SlackString(str::String)
    str2="{\"text\":\""*str*"\"}"
    run(`curl -X POST -H 'Content-type: application/json' --data $str2 $IWhU`;wait=false)
    return nothing
end

function SlackDict(dic::Dict)
    DIC=""
    for key ∈ sort(collect(keys(dic)))
        DIC=DIC*key*" : "*string(dic[key])*"\n"
    end
    DIC="{\"text\":\""*DIC*"\"}"
    run(`curl -X POST -H 'Content-type: application/json' --data $DIC $IWhU`;wait=false)
    return nothing
end

function SlackFile(filename;comment="")
    run(`curl -F file=@$filename -F "initial_comment=$comment" -F channels=$ChID -H "Authorization: Bearer $OAAT" https://slack.com/api/files.upload`;wait=false)
    return nothing
end
