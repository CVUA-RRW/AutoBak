rule flag:
    output:
        flag = "merge_flag",
    shell:
        "touch {output.flag}"