"""
    _pr_metafmt(level::Logging.LogLevel, _module, group, id, file, line)

MetaFormatter for ConsoleLogger for PR to adjust log message format
"""
function _pr_metafmt(level::Logging.LogLevel, _module, group, id, file, line)
    @nospecialize
    color = Logging.default_logcolor(level)
    prefix =
        "$(_module) | " *
        (level == Logging.Warn ? "Warning" : string(level)) *
        " ] :"
    suffix = ""
    Logging.Info <= level < Logging.Warn && return color, prefix, suffix
    _module !== nothing && (suffix *= "$(_module)")
    if file !== nothing
        _module !== nothing && (suffix *= " ")
        suffix *= Base.contractuser(file)
        if line !== nothing
            suffix *= ":$(isa(line, UnitRange) ? "$(first(line))-$(last(line))" : line)"
        end
    end
    !isempty(suffix) && (suffix = "@ " * suffix)

    return color, prefix, suffix
end

"""
    silence!()

Sets loglevel for PR to :Error, silencing Info and Warn
"""
function silence!()
    return set_logging_level!(:Error)
end

"""
    reset_logging_level!()

Resets the log level to Info
"""
function reset_logging_level!()
    Logging.global_logger(_LOGGER)

    return
end

"""
    restore_global_logger!()

Restores the global logger to its default state (before PR was loaded)
"""
function restore_global_logger!()
    Logging.global_logger(_DEFAULT_LOGGER)

    return
end

"""
    set_logging_level!(level::Symbol)

Sets the logging level for PR: :Info, :Warn, :Error, :Debug
"""
function set_logging_level!(level::Symbol)
    Logging.global_logger(_make_filtered_logger(getfield(Logging, level)))

    return
end

"""
    _make_filtered_logger(level::Logging.LogLevel)

Helper function to create the filtered logger for PR
"""
function _make_filtered_logger(level)
    LoggingExtras.EarlyFilteredLogger(_LOGGER) do log
        if log._module == PolyhedralRelaxations && log.level < level
            return false
        else
            return true
        end
    end
end
