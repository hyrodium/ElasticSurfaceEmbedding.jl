module Slack
    export SlackSettings, SlackString, SlackDict, SlackFile

    function SlackSettings(;IncomingWebhookURL="", OAuthAccessToken="", ChannelID="")
        global IWhU = IncomingWebhookURL
        global OAAT = OAuthAccessToken
        global ChID = ChannelID
        return nothing
    end

    function SlackString(str::String)
        str2="{\"text\":\""*str*"\"}"
        run(`curl -X POST -H 'Content-type: application/json' --data $str2 $IWhU`)
    end

    function SlackDict(dic::Dict)
        DIC=""
        for key in sort(collect(keys(dic)))
            DIC=DIC*key*" : "*string(dic[key])*"\n"
        end
        DIC="{\"text\":\""*DIC*"\"}"
        run(`curl -X POST -H 'Content-type: application/json' --data $DIC $IWhU`)
    end

    function SlackFile(filename;comment="")
        run(`curl -F file=@$filename -F "initial_comment=$comment" -F channels=$ChID -H "Authorization: Bearer $OAAT" https://slack.com/api/files.upload`)
    end
end
