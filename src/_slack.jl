struct SlackConfig
    channel::String
    token::String
end

SLACK = SlackConfig("","")

"""
    config_slack(;channel, token)

Set the channel and token for Slack bot. (optional)
"""
function config_slack(;channel, token)
    global SLACK = SlackConfig(channel, token)
end

function _send_file_to_slack(filename; comment = "")
    token = SLACK.token
    channel = SLACK.channel
    if token ≠ "" && channel ≠ ""
        if isfile(filename)
            cmd = `curl -F file=@$filename -F "initial_comment=$comment" -F channels=$(channel) -H "Authorization: Bearer $(token)" https://slack.com/api/files.upload`
        else
            cmd = `curl -d "text=$comment" -d "channel=$(channel)" -H "Authorization: Bearer $(token)" -X POST https://slack.com/api/chat.postMessage`
        end
        run(cmd, wait=false)
    end
    return
end
