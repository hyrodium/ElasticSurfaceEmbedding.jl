module Slack
    export SlackString, SlackDict, SlackFile

    function loadconfig()
        open("slack.txt") do file
            IWhU = readline(file)
            OAAT = readline(file)
            ChID = readline(file)
            return IWhU, OAAT, ChID
        end
    end

    function SlackString(str::String)
        IWhU, OAAT, ChID=loadconfig()
        str2="{\"text\":\""*str*"\"}"
        run(`curl -X POST -H 'Content-type: application/json' --data $str2 $IWhU`)
    end

    function SlackDict(dic::Dict)
        IWhU, OAAT, ChID=loadconfig()
        DIC=""
        for key ∈ sort(collect(keys(dic)))
            DIC=DIC*key*" : "*string(dic[key])*"\n"
        end
        DIC="{\"text\":\""*DIC*"\"}"
        run(`curl -X POST -H 'Content-type: application/json' --data $DIC $IWhU`)
    end

    function SlackFile(filename;comment="")
        IWhU, OAAT, ChID=loadconfig()
        run(`curl -F file=@$filename -F "initial_comment=$comment" -F channels=$ChID -H "Authorization: Bearer $OAAT" https://slack.com/api/files.upload`)
    end
end
