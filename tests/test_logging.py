"""Tests for the logging module (add_StreamHandler and add_FileHandler)."""

import logging
import tempfile
from pathlib import Path

import pytest

import metobs_toolkit


@pytest.fixture(autouse=True)
def cleanup_handlers():
    """Remove all handlers from the metobs_toolkit logger before and after each test."""
    logger = logging.getLogger("<metobs_toolkit>")
    logger.handlers.clear()
    yield
    logger.handlers.clear()


class TestAddStreamHandler:
    """Tests for add_StreamHandler function."""

    def test_add_stream_handler_default(self):
        """Test adding a StreamHandler with default settings."""
        metobs_toolkit.add_StreamHandler()

        logger = logging.getLogger("<metobs_toolkit>")
        stream_handlers = [
            h
            for h in logger.handlers
            if isinstance(h, logging.StreamHandler)
            and not isinstance(h, logging.FileHandler)
        ]

        assert len(stream_handlers) == 1
        assert stream_handlers[0].level == logging.DEBUG

    def test_add_stream_handler_custom_level(self):
        """Test adding a StreamHandler with custom log level."""
        metobs_toolkit.add_StreamHandler(setlvl="WARNING")

        logger = logging.getLogger("<metobs_toolkit>")
        stream_handlers = [
            h
            for h in logger.handlers
            if isinstance(h, logging.StreamHandler)
            and not isinstance(h, logging.FileHandler)
        ]

        assert len(stream_handlers) == 1
        assert stream_handlers[0].level == logging.WARNING

    def test_no_duplicate_stream_handler(self):
        """Test that duplicate StreamHandlers are not added."""
        metobs_toolkit.add_StreamHandler(setlvl="DEBUG")
        metobs_toolkit.add_StreamHandler(setlvl="INFO")

        logger = logging.getLogger("<metobs_toolkit>")
        stream_handlers = [
            h
            for h in logger.handlers
            if isinstance(h, logging.StreamHandler)
            and not isinstance(h, logging.FileHandler)
        ]

        # Should only have one handler (first one at DEBUG covers INFO)
        assert len(stream_handlers) == 1


class TestAddFileHandler:
    """Tests for add_FileHandler function."""

    def test_add_file_handler_creates_file(self):
        """Test that add_FileHandler creates the log file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            logfile = Path(tmpdir) / "test.log"

            metobs_toolkit.add_FileHandler(filepath=logfile)

            assert logfile.exists()

    def test_add_file_handler_writes_logs(self):
        """Test that FileHandler writes logs to the file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            logfile = Path(tmpdir) / "test.log"

            metobs_toolkit.add_FileHandler(filepath=logfile, setlvl="INFO")

            logger = logging.getLogger("<metobs_toolkit>")
            logger.info("Test log message")

            # Flush handlers
            for handler in logger.handlers:
                handler.flush()

            content = logfile.read_text()
            assert "Test log message" in content

    def test_add_file_handler_custom_level(self):
        """Test adding a FileHandler with custom log level."""
        with tempfile.TemporaryDirectory() as tmpdir:
            logfile = Path(tmpdir) / "test.log"

            metobs_toolkit.add_FileHandler(filepath=logfile, setlvl="ERROR")

            logger = logging.getLogger("<metobs_toolkit>")
            file_handlers = [
                h for h in logger.handlers if isinstance(h, logging.FileHandler)
            ]

            assert len(file_handlers) == 1
            assert file_handlers[0].level == logging.ERROR

    def test_clearlog_overwrites_file(self):
        """Test that clearlog=True clears the log file."""
        with tempfile.TemporaryDirectory() as tmpdir:
            logfile = Path(tmpdir) / "test.log"

            # Create file with existing content
            logfile.write_text("Existing content\n")

            metobs_toolkit.add_FileHandler(filepath=logfile, clearlog=True)

            content = logfile.read_text()
            assert "Existing content" not in content
