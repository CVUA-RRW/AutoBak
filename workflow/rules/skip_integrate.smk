rule flag:
    output:
        flag = "flags/merge_flag",
    shell:
        "touch {output.flag}"