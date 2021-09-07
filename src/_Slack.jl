function _send_file_to_slack(filename; comment = "")
    cmd = `curl -F file=@$filename -F "initial_comment=$comment" -F channels=$ChID -H "Authorization: Bearer $OAAT" https://slack.com/api/files.upload`
    run(cmd; wait = false)
    return
end
