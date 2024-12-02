process PUBLISH_COMMAND {

    /* */

    publishDir params.run_command, mode: "copy", overwrite: true

    errorStrategy { task.attempt < 2 ? 'retry' : 'ignore' }
    maxRetries 1

    cpus 1

    output:
    path "run_command.sh"

    script:
    command = workflow.commandLine
    """
    echo "${command}" > run_command.sh
    """
}
