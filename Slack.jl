module Slack
    export SlackSetting, SlackString, SlackDict, SlackFile

    function SlackSetting(val)
        global setting=val
        return nothing
    end

    function SlackString(str::String)
        IncomingWebhookURL=setting["IncomingWebhookURL"]
        str2="{\"text\":\""*str*"\"}"
        run(`curl -X POST -H 'Content-type: application/json' --data $str2 $IncomingWebhookURL`)
    end

    function SlackDict(dic::Dict)
        DIC=""
        for key in sort(collect(keys(dic)))
            DIC=DIC*key*" : "*string(dic[key])*"\n"
        end
        DIC="{\"text\":\""*DIC*"\"}"
        IncomingWebhookURL=setting["IncomingWebhookURL"]
        run(`curl -X POST -H 'Content-type: application/json' --data $DIC $IncomingWebhookURL`)
    end

    function SlackFile(filename;comment="")
        ChannelID=setting["ChannelID"]
        OAuthAccessToken=setting["OAuthAccessToken"]
        run(`curl -F file=@$filename -F "initial_comment=$comment" -F channels=$ChannelID -H "Authorization: Bearer $OAuthAccessToken" https://slack.com/api/files.upload`)
    end
end
